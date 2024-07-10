#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include "pair_pwmlff.h"
#include "nep_cpu.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"
#include "domain.h"
#include <dlfcn.h>
#include <cuda_runtime.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPWMLFF::PairPWMLFF(LAMMPS *lmp) : Pair(lmp)
{
    me = comm->me;
	writedata = 1;
    comm_reverse = 3;

}

PairPWMLFF::~PairPWMLFF()
{
    if (allocated) 
    {
        memory->destroy(setflag);
        memory->destroy(cutsq);
        memory->destroy(f_n);
        memory->destroy(e_atom_n);
    }

    if (me == 0) {
        fclose(explrError_fp);
    }

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPWMLFF::allocate()
{
    allocated = 1;
    int np1 = atom->ntypes ;
    memory->create(setflag, np1 + 1, np1 + 1, "pair:setflag");
    for (int i = 1; i <= np1; i++)
        for (int j = i; j <= np1; j++) setflag[i][j] = 0;
    memory->create(cutsq, np1 + 1, np1 + 1, "pair:cutsq");

}

/* ----------------------------------------------------------------------
   global settings pair_style 
------------------------------------------------------------------------- */

void PairPWMLFF::settings(int narg, char** arg)
{
    int ff_idx;
    int iarg = 1;  // index of arg after 'num_ff'
    int rank;
    int num_devices;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    num_devices = torch::cuda::device_count();

    if (narg <= 0) error->all(FLERR, "Illegal pair_style command"); // numbers of args after 'pair_style pwmlff'
    std::vector<std::string> models;

    num_ff = utils::inumeric(FLERR, arg[0], false, lmp);    // number of models
    // for (int ii = 1; ii < iarg; ++ii) {
    //     models.push_back(arg[ii]);                          // model files
    // }
    for (int ii = 0; ii < num_ff; ++ii) {
        models.push_back(arg[iarg++]);                          // model files
    }
    while (iarg < narg) {
        if (strcmp(arg[iarg], "out_freq") == 0) {
            out_freq = utils::inumeric(FLERR, arg[++iarg], false, lmp);
        } else if (strcmp(arg[iarg], "out_file") == 0) {
            explrError_fname = arg[++iarg];
        } 
        iarg++;
    }

    if (me == 0) {
        explrError_fp = fopen(&explrError_fname[0], "w");
    }

    // device = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;
    torch::DeviceType device_type = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;
    device = device_type == torch::kCUDA ?  torch::Device(device_type, rank % num_devices) : torch::Device(device_type);
    dtype = torch::kFloat64;

    std::cout<<"the numbor of gpu is " <<  num_devices << std::endl;
    std::cout<<"the mpi process is " << comm->nprocs << std::endl;
    
    if (me == 0) utils::logmesg(this -> lmp, "<---- Loading model ---->");
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        std::string model_file = models[ff_idx];
        try
        {   // load jit model, they are DP, or NEP
            // module = torch::jit::load(model_file, c10::Device(device)); 
            module = torch::jit::load(model_file, device);
            module.to(dtype);
            // module.eval();
            modules.push_back(module);
            if (me == 0) printf("\nLoading jitscript model file:   %s\n", model_file.c_str());
            model_type = 0;
        }
        catch (const c10::Error e)
        {
            // load txt format model file, it should be nep
            try {
                bool is_rank_0 = (comm->me == 0);
                void* handle = dlopen("libnep_gpu.so", RTLD_LAZY);
                if (!handle) {
                std::cout << "Could not find libnep_gpu.so " << dlerror() << " The nep cpu version will be used for lammps! " << std::endl;
                use_nep_gpu = false;
                } else if (num_devices < 1) {
                    std::cout << "Can not find the GPU devices. The nep cpu version will be used for lammps! " << std::endl;
                    use_nep_gpu = false;
                    dlclose(handle);
                }
                else {
                    std::cout << "Load libnep_gpu.so success. The GPU device will be used !"<< std::endl;
                    use_nep_gpu = true;
                    dlclose(handle);
                }
                // if the gpu nums > 0 and libnep.so is exsits, use gpu model
                if (use_nep_gpu) {
                    int device_id = rank % num_devices;
                    cudaSetDevice(device_id);
                    nep_gpu_model.init_from_file(model_file.c_str(), is_rank_0, device_id);
                    model_type = 2;
                    printf("MPI rank %d rank using GPU device %d\n", rank, device_id);
                    // std::cout<<"load nep.txt success and the model type is 2" << std::endl;
                } else {
                    nep_cpu_model.init_from_file(model_file, is_rank_0);
                    nep_cpu_models.push_back(nep_cpu_model);
                    model_type = 1;
                }
                if (me == 0) printf("\nLoading txt model file:   %s\n", model_file.c_str());
                
            } catch (const c10::Error e) {
                std::cerr << "Failed to load model :" << e.msg() << std::endl;
            }
        }
    }
    if (model_type == 0) {
        if ((device_type == torch::kCUDA) && (rank == 0)) {
            if (num_devices < comm->nprocs) {
            std::cout << "----------------------------------------------------------------------------------" << std::endl;
            std::cout << " Warning: There are " << num_devices << " GPUs available " << std::endl;
            std::cout << " But have " << comm->nprocs << " MPI processes, may result in poor performance!!!" << std::endl;
            std::cout << "----------------------------------------------------------------------------------" << std::endl;
            }
        }
        try {
            torch::jit::IValue model_type_value = module.attr("model_type");
            model_name = model_type_value.toString()->string();
        }
        catch (const torch::jit::ObjectAttributeError& e) {
            model_name = "DP";
        }
        if (model_name == "DP") {
            cutoff = module.attr("Rmax").toDouble();
        } else if (model_name == "NEP") {
            // 设置nep的参数
            cutoff_radial  = module.attr("cutoff_radial").toDouble();
            cutoff = cutoff_radial;
            cutoff_angular = module.attr("cutoff_angular").toDouble();
        } else {
            std::cout << "ERROR: the model_type of input model " << model_name << " is not supported! Please check the input model! " << std::endl;
            error->universe_all(FLERR, "ERROR: the model_type of input model is not supported! Please check the input model!");
        }

        //common params of DP and NEP
        max_neighbor = module.attr("maxNeighborNum").toInt();
        
        // print information
        if (me == 0) {
        utils::logmesg(this -> lmp, "<---- Load model successful!!! ---->");
        printf("\nDevice:       %s", device == torch::kCPU ? "CPU" : "GPU");
        printf("\nModel type:   %5d",5);
        printf("\nModel nums:   %5d",num_ff);
        printf("\ncutoff :      %12.6f",cutoff);
        printf("\nmax_neighbor: %5d\n", max_neighbor);
        }
    } else if (model_type == 1) {
        cutoff = nep_cpu_model.paramb.rc_radial;
    } else if (model_type == 2) {
        cutoff = nep_gpu_model.paramb.rc_radial;
    }
    // since we need num_ff, so well allocate memory here
    // but not in allocate()
    memory->create(f_n, num_ff, atom->nmax, 3, "pair_pwmlff:f_n");
    memory->create(e_atom_n, num_ff, atom->natoms, "pair_pwmlff:e_atom_n");
} 

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs pair_coeff 
------------------------------------------------------------------------- */

void PairPWMLFF::coeff(int narg, char** arg)
{
    int ntypes = atom->ntypes;
    if (!allocated) { allocate(); }

    // pair_coeff * * 
    int ilo, ihi, jlo, jhi;
    utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error); // arg[0] = *
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error); // arg[1] = *

    int count = 0;
    for(int i = ilo; i <= ihi; i++) {
        for(int j = MAX(jlo,i); j <= jhi; j++) 
        {
            setflag[i][j] = 1;
            count++;
        }
    }

    if (model_type == 0) {
        auto atom_type_module = module.attr("atom_type").toList();
        model_ntypes = atom_type_module.size();
        if (ntypes > model_ntypes || ntypes != narg - 2)  // type numbers in strucutre file and in pair_coeff should be the same
        {
            error->all(FLERR, "Element mapping is not correct, ntypes = " + std::to_string(ntypes));
        }
        for (int ii = 2; ii < narg; ++ii) {
            int temp = std::stoi(arg[ii]);
            auto iter = std::find(atom_type_module.begin(), atom_type_module.end(), temp);   
            if (iter != atom_type_module.end() || arg[ii] == 0)
            {
                int index = std::distance(atom_type_module.begin(), iter);
                model_atom_type_idx.push_back(index);
                atom_types.push_back(temp);
            }
            else
            {
                error->all(FLERR, "This element is not included in the machine learning force field");
            }
        }
    } else if (model_type == 1) {
        std::vector<int> atom_type_module = nep_cpu_model.element_atomic_number_list;
        model_ntypes = atom_type_module.size();
        if (ntypes > model_ntypes || ntypes != narg - 2)  // type numbers in strucutre file and in pair_coeff should be the same
        {
            error->all(FLERR, "Element mapping is not correct, ntypes = " + std::to_string(ntypes));
        }
        for (int ii = 2; ii < narg; ++ii) {
            int temp = std::stoi(arg[ii]);
            auto iter = std::find(atom_type_module.begin(), atom_type_module.end(), temp);   
            if (iter != atom_type_module.end() || arg[ii] == 0)
            {
                int index = std::distance(atom_type_module.begin(), iter);
                model_atom_type_idx.push_back(index); 
                atom_types.push_back(temp);
                for(int jj=0; jj < num_ff; ++jj){
                    nep_cpu_models[jj].map_atom_types.push_back(temp);
                    nep_cpu_models[jj].map_atom_type_idx.push_back(index);
                }
                // std::cout<<"=== the config atom type "<< temp << " index in ff is "  << index << std::endl;
            }
            else
            {
                error->all(FLERR, "This element is not included in the machine learning force field");
            }
        }
    } else if (model_type == 2) { // for nep_gpu
        // check or reset
        std::vector<int> atom_type_module = nep_gpu_model.element_atomic_number_list;
        model_ntypes = atom_type_module.size();
        if (ntypes > model_ntypes || ntypes != narg - 2)  // type numbers in strucutre file and in pair_coeff should be the same
        {
            error->all(FLERR, "Element mapping is not correct, ntypes = " + std::to_string(ntypes));
        }
        for (int ii = 2; ii < narg; ++ii) {
            int temp = std::stoi(arg[ii]);
            auto iter = std::find(atom_type_module.begin(), atom_type_module.end(), temp);   
            if (iter != atom_type_module.end() || arg[ii] == 0)
            {
                int index = std::distance(atom_type_module.begin(), iter);
                model_atom_type_idx.push_back(index); 
                atom_types.push_back(temp);
                    nep_gpu_model.map_atom_type_idx.push_back(index);
                // std::cout<<"=== the config atom type "<< temp << " index in ff is "  << index << std::endl;
            }
            else
            {
                error->all(FLERR, "This element is not included in the machine learning force field");
            }
        }        
    }
   if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairPWMLFF::init_one(int i, int j)
{
    //if (setflag[i][j] == 0) { error->all(FLERR, "All pair coeffs are not set"); 

    return cutoff;
}


void PairPWMLFF::init_style()
{
    if (force->newton_pair == 0) error->all(FLERR, "Pair style PWMLFF requires newton pair on");
    // Using a nearest neighbor table of type full
    neighbor->add_request(this, NeighConst::REQ_FULL);
}
/* ---------------------------------------------------------------------- */

std::pair<double, double> PairPWMLFF::calc_max_error(double ***f_n, double **e_atom_n)
{
    int i, j;
    int ff_idx;
    double max_err, err, max_err_ei, err_ei;
    double num_ff_inv;
    int nlocal = atom->nlocal;
    // int *tag = atom->tag;

    max_err = -1.0;
    max_err_ei = -1.0;
    num_ff_inv = 1.0 / num_ff;

    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        p_ff_idx = ff_idx;
        comm->reverse_comm(this);
    }

    std::vector<double> f_ave;
    std::vector<double> f_err[num_ff];
    std::vector<double> ei_ave;
    std::vector<double> ei_err[num_ff];

    f_ave.resize(nlocal * 3);
    ei_ave.resize(nlocal);

    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        f_err[ff_idx].resize(nlocal * 3);
        ei_err[ff_idx].resize(nlocal);
    }

    // sum over all models
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        for (i = 0; i < nlocal; i++) {
            // std::cout << "f_n[" << ff_idx << "][" << i << "][0] = " << tag[i] << " " << f_n[ff_idx][i][0] << std::endl;
            f_ave[i * 3 + 0] += f_n[ff_idx][i][0];
            f_ave[i * 3 + 1] += f_n[ff_idx][i][1];
            f_ave[i * 3 + 2] += f_n[ff_idx][i][2];
            ei_ave[i] += e_atom_n[ff_idx][i];
            // std::cout<< "ff " << ff_idx << " i " << i << " ei " <<  e_atom_n[ff_idx][i] << " force " << f_n[ff_idx][i][0] << " "  << f_n[ff_idx][i][1] << " "  << f_n[ff_idx][i][2] << std::endl;
        }
    }

    // calc ensemble average
    for (i = 0; i < 3 * nlocal; i++) {
        f_ave[i] *= num_ff_inv;
    }
    for (i = 0; i < nlocal; i++) {
        ei_ave[i] *= num_ff_inv;
    }

    // calc error
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        for (i = 0; i < nlocal; i++) {
            f_err[ff_idx][i * 3 + 0] = f_n[ff_idx][i][0] - f_ave[i * 3 + 0];
            f_err[ff_idx][i * 3 + 1] = f_n[ff_idx][i][1] - f_ave[i * 3 + 1];
            f_err[ff_idx][i * 3 + 2] = f_n[ff_idx][i][2] - f_ave[i * 3 + 2];
            ei_err[ff_idx][i] = e_atom_n[ff_idx][i] - ei_ave[i];
        }
    }

    // find max error
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        for (j = 0; j < nlocal * 3; j += 3) {
            err = f_err[ff_idx][j] * f_err[ff_idx][j] + f_err[ff_idx][j + 1] * f_err[ff_idx][j + 1] + f_err[ff_idx][j + 2] * f_err[ff_idx][j + 2];
            err = sqrt(err);
            if (err > max_err) max_err = err;
        }
        for (j = 0; j < nlocal; j++) {
            err_ei = ei_err[ff_idx][j];
            if (err_ei > max_err_ei) max_err_ei = err_ei;
        }
    }
    return std::make_pair(max_err, max_err_ei);

}

int PairPWMLFF::pack_reverse_comm(int n, int first, double* buf) {
    int i, m, last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
        buf[m++] = f_n[p_ff_idx][i][0];
        buf[m++] = f_n[p_ff_idx][i][1];
        buf[m++] = f_n[p_ff_idx][i][2];
    }
    return m;
}

void PairPWMLFF::unpack_reverse_comm(int n, int* list, double* buf) {
    int i, j, m;

    m = 0;
    for (i = 0; i < n; i++) {
        j = list[i];
        f_n[p_ff_idx][j][0] += buf[m++];
        f_n[p_ff_idx][j][1] += buf[m++];
        f_n[p_ff_idx][j][2] += buf[m++];
    }

}

void PairPWMLFF::grow_memory()
{
  if (atom->nmax > nmax) {
    printf("@@@ allocate new %7d %7d %7d\n", update->ntimestep, nmax, atom->nmax);
    nmax = atom->nmax;
    memory->grow(f_n, num_ff, nmax, 3, "pair_pwmlff:f_n");
    memory->grow(e_atom_n, num_ff, nmax, "pair_pwmlff:e_atom_n");
  }
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<double>> PairPWMLFF::generate_neighdata()
{   
    int i, j, k, ii, jj, inum, jnum, itype, jtype;
    double xtmp, ytmp, ztmp, delx, dely, delz, rsq, rij;
    int *ilist, *jlist, *numneigh, **firstneigh;
    int etnum;

    double **x = atom->x;
    int *type = atom->type;
    int *tag = atom->tag;
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    // int ntypes = atom->ntypes;
    int ntypes = model_ntypes;
    int n_all = nlocal + nghost;
    double rc2 = cutoff * cutoff;

    double min_dR = 1000;
    double min_dR_all;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    std::vector<std::vector<int>> num_neigh(inum, std::vector<int>(ntypes));
    // imagetype.resize(inum);
    imagetype_map.resize(inum);
    neighbor_list.resize(inum * ntypes * max_neighbor);
    dR_neigh.resize(inum * ntypes * max_neighbor * 4);
    // use_type.resize(n_all);
    std::vector<int> type_to_model(n_all);


    for (int ii = 0; ii < n_all; ii++)
    {
        // use_type[ii] = atom_types[type[ii] - 1];
        type_to_model[ii] = model_atom_type_idx[type[ii] - 1] + 1;
        // type[0], type[1], type[2], type[3], type[4], : 2, 2, 1, 2, 2, ...
        // atom_types[0], atom_types[1] : 6, 1
        // use_type[0], use_type[1], use_type[2], use_type[3], use_type[4] : 1, 1, 6, 1, 1, ...
    }
    for (i = 0; i < nlocal; i++) {
        for (j = 0; j < ntypes; j++) {
            num_neigh[i][j] = 0;
        }
    }
    for (ii = 0; ii < inum; ii++)               // local atoms: 5, CH4
    {    
        i = ilist[ii];                          // 0, 1, 2, 3, 4
        // itype = type[i];                        // 2, 2, 1, 2, 2
        itype = type_to_model[i];                   // 1, 1, 3, 1, 1
        jlist = firstneigh[i];
        jnum = numneigh[i];                     // 4, 4, 4, 4, 4
        imagetype_map[ii] = itype - 1;          // 1, 1, 0, 1, 1        python index from 0
        // imagetype[ii] = use_type[i];            // 1, 1, 6, 1, 1

        for (jj = 0; jj < jnum; jj++)
        {
            j = jlist[jj];                      // 1, 2, 3, 4;   0, 2, 3, 4;   0, 1, 3, 4;   0, 1, 2, 4;   0, 1, 2, 3
            delx = x[j][0] - x[i][0];
            dely = x[j][1] - x[i][1];
            delz = x[j][2] - x[i][2];
            rsq = delx * delx + dely * dely + delz * delz;
            // jtype = type[j];                    // 2, 1, 2, 2;   2, 1, 2, 2;   2, 2, 2, 2;   2, 2, 1, 2;   2, 2, 1, 2
            jtype = type_to_model[j];               // 1, 3, 1, 1;   1, 3, 1, 1;   1, 1, 1, 1;   1, 1, 3, 1;   1, 1, 3, 1
            if (rsq <= rc2) 
            {
                etnum = num_neigh[i][jtype - 1];
                rij = sqrt(rsq);
                int index = i * ntypes * max_neighbor + (jtype - 1) * max_neighbor + etnum;
                dR_neigh[index * 4 + 0] = rij;
                dR_neigh[index * 4 + 1] = delx;
                dR_neigh[index * 4 + 2] = dely;
                dR_neigh[index * 4 + 3] = delz;
                neighbor_list[index] = j + 1;
                num_neigh[i][jtype - 1] += 1;
                // std::cout << "num_neigh[" << i << "][" << jtype - 1 << "] = " << num_neigh[i][jtype - 1] << std::endl;
                if (rsq < min_dR) min_dR = rsq;
                // std::cout<< "nlist index " << index << " j " << j+1 << std::endl;
            }
        }
    }

    MPI_Allreduce(&min_dR, &min_dR_all, 1, MPI_DOUBLE, MPI_MIN, world);

    if (min_dR_all < 0.81) {
        if (me == 0) {
            std::cout << "ERROR: there are two atoms too close, min_dR_all = " << min_dR_all << std::endl;
        }
        error->universe_all(FLERR, "there are two atoms too close");
    }
    return std::make_tuple(std::move(imagetype_map), std::move(neighbor_list), std::move(dR_neigh));
    // return std::make_tuple(imagetype, imagetype_map, neighbor_list, dR_neigh);
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<double>> PairPWMLFF::generate_neighdata_nep()
{   
    int i, j, k, ii, jj, inum, jnum, itype, jtype;
    double xtmp, ytmp, ztmp, delx, dely, delz, rsq, rij;
    int *ilist, *jlist, *numneigh, **firstneigh;
    int etnum;

    double **x = atom->x;
    int *type = atom->type;
    int *tag = atom->tag;
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    // int ntypes = atom->ntypes;
    int ntypes = model_ntypes;
    int n_all = nlocal + nghost;
    double rc2 = cutoff * cutoff;

    double min_dR = 1000;
    double min_dR_all;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    std::vector<std::vector<int>> num_neigh(inum, std::vector<int>(ntypes));
    // imagetype.resize(inum);
    imagetype_map.assign(inum, -1);
    neighbor_list.assign(inum * ntypes * max_neighbor, 0);
    neighbor_type_list.assign(inum * ntypes * max_neighbor, -1);
    dR_neigh.assign(inum * ntypes * max_neighbor * 4, 0);
    // use_type.resize(n_all);
    std::vector<int> type_to_model(n_all);


    for (int ii = 0; ii < n_all; ii++)
    {
        // use_type[ii] = atom_types[type[ii] - 1];
        type_to_model[ii] = model_atom_type_idx[type[ii] - 1] + 1;
        // type[0], type[1], type[2], type[3], type[4], : 2, 2, 1, 2, 2, ...
        // atom_types[0], atom_types[1] : 6, 1
        // use_type[0], use_type[1], use_type[2], use_type[3], use_type[4] : 1, 1, 6, 1, 1, ...
    }
    for (i = 0; i < nlocal; i++) {
        for (j = 0; j < ntypes; j++) {
            num_neigh[i][j] = 0;
        }
    }
    for (ii = 0; ii < inum; ii++)               // local atoms: 5, CH4
    {    
        i = ilist[ii];                          // 0, 1, 2, 3, 4
        // itype = type[i];                        // 2, 2, 1, 2, 2
        itype = type_to_model[i];                   // 1, 1, 3, 1, 1
        jlist = firstneigh[i];
        jnum = numneigh[i];                     // 4, 4, 4, 4, 4
        imagetype_map[ii] = itype - 1;          // 1, 1, 0, 1, 1        python index from 0
        // imagetype[ii] = use_type[i];            // 1, 1, 6, 1, 1
        int count_radial = 0;
        for (jj = 0; jj < jnum; jj++)
        {
            j = jlist[jj];                      // 1, 2, 3, 4;   0, 2, 3, 4;   0, 1, 3, 4;   0, 1, 2, 4;   0, 1, 2, 3
            delx = x[j][0] - x[i][0];
            dely = x[j][1] - x[i][1];
            delz = x[j][2] - x[i][2];
            rsq = delx * delx + dely * dely + delz * delz;
            // jtype = type[j];                    // 2, 1, 2, 2;   2, 1, 2, 2;   2, 2, 2, 2;   2, 2, 1, 2;   2, 2, 1, 2
            jtype = type_to_model[j];               // 1, 3, 1, 1;   1, 3, 1, 1;   1, 1, 1, 1;   1, 1, 3, 1;   1, 1, 3, 1
            if (rsq <= rc2) 
            {   
                //Wuxing: the DP use i to set dR_neigh, but I think it should be ii, the output of ilist is 0, 1, 2, 3, 4... , same as i 

                // etnum = num_neigh[i][jtype - 1]; 
                // rij = sqrt(rsq);
                // int index = i * ntypes * max_neighbor + (jtype - 1) * max_neighbor + etnum;
                // dR_neigh[index * 4 + 0] = rij;
                // dR_neigh[index * 4 + 1] = delx;
                // dR_neigh[index * 4 + 2] = dely;
                // dR_neigh[index * 4 + 3] = delz;
                // neighbor_list[index] = j + 1;
                // num_neigh[i][jtype - 1] += 1;
                // // std::cout << "num_neigh[" << i << "][" << jtype - 1 << "] = " << num_neigh[i][jtype - 1] << std::endl;
                // if (rsq < min_dR) min_dR = rsq;

                rij = sqrt(rsq);
                int index = i * ntypes * max_neighbor;
                neighbor_list[index + count_radial ] = j + 1;
                neighbor_type_list[index + count_radial] = jtype-1;
                dR_neigh[index*4 + 4*count_radial  ] = rij;
                dR_neigh[index*4 + 4*count_radial+1] = delx;
                dR_neigh[index*4 + 4*count_radial+2] = dely;
                dR_neigh[index*4 + 4*count_radial+3] = delz;
                count_radial++;
                // std::cout<< "atom_I_ii " << ii << " index_I " << i << " type_I " << itype << "    atom_j_jj " << jj << " index_j " << j << " type_j " << jtype << " neigbor_start "<< index << std::endl;

                // std::cout<< "nlist index " << index + count_radial << " count " << count_radial << " j " << j+1 << std::endl;

            }
        }
    }

    MPI_Allreduce(&min_dR, &min_dR_all, 1, MPI_DOUBLE, MPI_MIN, world);

    if (min_dR_all < 0.81) {
        if (me == 0) {
            std::cout << "ERROR: there are two atoms too close, min_dR_all = " << min_dR_all << std::endl;
        }
        error->universe_all(FLERR, "there are two atoms too close");
    }
    return std::make_tuple(std::move(imagetype_map), std::move(neighbor_list), std::move(neighbor_type_list), std::move(dR_neigh));
    // return std::make_tuple(imagetype, imagetype_map, neighbor_list, dR_neigh);
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<double>>PairPWMLFF::convert_dim(bool is_build_neighbor){
    int nlocal = atom->nlocal;
    int nghost = atom->nghost; 
    int n_all = nlocal + nghost;
    int *itype, *numneigh, **firstneigh;
    int ii, jj, inum, jnum;
    
    inum = list->inum;
    itype = atom->type;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    double **x = atom->x;
    itype_convert_map.assign(n_all, -1);
    position_cpu.assign(n_all*3, 0);
    firstneighbor_cpu.assign(inum * nep_gpu_nm, 0);

    for(ii=0; ii < n_all;ii++) {
        itype_convert_map[ii] = model_atom_type_idx[itype[ii] - 1];
        position_cpu[ii] = x[ii][0];
        position_cpu[  n_all+ii] = x[ii][1];
        position_cpu[2*n_all+ii] = x[ii][2];
        if (ii < inum){
            jnum = numneigh[ii];
            for(jj=0; jj < jnum; ++jj){
                firstneighbor_cpu[ii*nep_gpu_nm + jj] = firstneigh[ii][jj];
            }
        }
    }
    return std::make_tuple(std::move(itype_convert_map), std::move(firstneighbor_cpu), std::move(position_cpu));
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<float>> PairPWMLFF::generate_neighdata_nep_gpu()
{   
    int i, j, k, ii, jj, inum, jnum, itype, jtype;
    double xtmp, ytmp, ztmp, delx, dely, delz, rsq, rij;
    int *ilist, *jlist, *numneigh, **firstneigh;
    int etnum;

    double **x = atom->x;
    int *type = atom->type;
    int *tag = atom->tag;
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    // int ntypes = atom->ntypes;
    int ntypes = model_ntypes;
    int n_all = nlocal + nghost;
    double rc2 = nep_gpu_model.paramb.rc_radial*nep_gpu_model.paramb.rc_radial;
    double rc_a2 = nep_gpu_model.paramb.rc_angular*nep_gpu_model.paramb.rc_angular;
    // std::cout << "========= rc2 " << rc2 << " rc_a2 " << rc_a2 << std::endl;
    double min_dR = 1000;
    double min_dR_all;

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    int NM = nep_gpu_nm;
    int N_NM = inum * NM;
    // imagetype.resize(inum);
    // printf("find neighbor n_all %d nlocal %d nghost %d inum %d\n", n_all, nlocal, nghost, inum);
    neighbor_list.assign(N_NM, 0);
    neighbor_angular_list.assign(N_NM, 0);
    itype_convert_map.assign(n_all, -1);
    neigbor_num_list.assign(inum, 0);
    neigbor_angular_num_list.assign(inum, 0);
    rij_nep_gpu.assign(N_NM * 6, 0);
    // use_type.resize(n_all);
    // std::vector<int> type_to_model(n_all);

    for (int ii = 0; ii < n_all; ii++)
    {
        itype_convert_map[ii] = model_atom_type_idx[type[ii] - 1];
    }

    for (ii = 0; ii < inum; ii++)               // local atoms: 5, CH4
    {    
        i = ilist[ii];                          // 0, 1, 2, 3, 4, Although we found i==ilist [ii] during testing, it is possible that there may not be enough testing
        // itype = type[i];                        // 2, 2, 1, 2, 2
        itype = itype_convert_map[i];                   // 1, 1, 3, 1, 1
        jlist = firstneigh[i];
        jnum = numneigh[i];                     // 4, 4, 4, 4, 4
        int count_radial = 0;
        int count_angular =0;
        for (jj = 0; jj < jnum; jj++)
        {
            j = jlist[jj];                      // 1, 2, 3, 4;   0, 2, 3, 4;   0, 1, 3, 4;   0, 1, 2, 4;   0, 1, 2, 3
            delx = x[j][0] - x[i][0];
            dely = x[j][1] - x[i][1];
            delz = x[j][2] - x[i][2];
            rsq = delx * delx + dely * dely + delz * delz;
            // jtype = type[j];                    // 2, 1, 2, 2;   2, 1, 2, 2;   2, 2, 2, 2;   2, 2, 1, 2;   2, 2, 1, 2
            jtype = itype_convert_map[j];               // 1, 3, 1, 1;   1, 3, 1, 1;   1, 1, 1, 1;   1, 1, 3, 1;   1, 1, 3, 1
            if (rsq <= rc2) 
            {
                neighbor_list[count_radial * inum + i] = j;
                rij_nep_gpu[0* N_NM + count_radial*inum+i] = delx;
                rij_nep_gpu[1* N_NM + count_radial*inum+i] = dely;
                rij_nep_gpu[2* N_NM + count_radial*inum+i] = delz;
                count_radial++;
            }
            if (rsq <= rc_a2) 
            {   
                neighbor_angular_list[count_angular * inum + i] = j;
                rij_nep_gpu[3* N_NM + count_angular*inum+i] = delx;
                rij_nep_gpu[4* N_NM + count_angular*inum+i] = dely;
                rij_nep_gpu[5* N_NM + count_angular*inum+i] = delz;
                count_angular++;
            }
        }
        neigbor_num_list[i] = count_radial;
        neigbor_angular_num_list[i] = count_angular;
    }

    MPI_Allreduce(&min_dR, &min_dR_all, 1, MPI_DOUBLE, MPI_MIN, world);

    if (min_dR_all < 0.81) {
        if (me == 0) {
            std::cout << "ERROR: there are two atoms too close, min_dR_all = " << min_dR_all << std::endl;
        }
        error->universe_all(FLERR, "there are two atoms too close");
    }
    return std::make_tuple(std::move(itype_convert_map), std::move(neighbor_list), std::move(neighbor_angular_list), std::move(neigbor_num_list), std::move(neigbor_angular_num_list), std::move(rij_nep_gpu));
    // return std::make_tuple(imagetype, imagetype_map, neighbor_list, dR_neigh);
}

void PairPWMLFF::compute(int eflag, int vflag)
{
    if (eflag || vflag) ev_setup(eflag, vflag);

    // int newton_pair = force->newton_pair;
    int ff_idx;
    int nlocal = atom->nlocal;
    int current_timestep = update->ntimestep;
    // int total_timestep = update->laststep;
    bool calc_virial_from_mlff = false;
    bool calc_egroup_from_mlff = false;
    int ntypes = atom->ntypes;
    int nghost = atom->nghost;
    int n_all = nlocal + nghost;
    // int inum, jnum, itype, jtype;
    double max_err, global_max_err, max_err_ei, global_max_err_ei;
    bool is_build_neighbor = (current_timestep % neighbor->every == 0);
    // for dp and nep model from jitscript
    if (model_type == 0) {

        int inum = list->inum;
        double *virial = force->pair->virial;
        double **x = atom->x;
        double **f = atom->f;
        int *type = atom->type;
        // auto t4 = std::chrono::high_resolution_clock::now();
        std::vector<int> imagetype_map, neighbor_list, neighbor_type_list;
        std::vector<double> dR_neigh;

        if (model_name == "DP") {
            std::tie(imagetype_map, neighbor_list, dR_neigh) = generate_neighdata();
        } else if (model_name == "NEP") {
            // std::cout<< "==== do nep find neigh ====" << std::endl;
            std::tie(imagetype_map, neighbor_list, neighbor_type_list, dR_neigh) = generate_neighdata_nep();
        }
        if (inum == 0) return;
        // auto t5 = std::chrono::high_resolution_clock::now();
        auto int_tensor_options = torch::TensorOptions().dtype(torch::kInt);
        auto float_tensor_options = torch::TensorOptions().dtype(torch::kFloat64);
        torch::Tensor imagetype_map_tensor = torch::from_blob(imagetype_map.data(), {inum}, int_tensor_options).to(device);
        torch::Tensor neighbor_list_tensor = torch::from_blob(neighbor_list.data(), {1, inum, max_neighbor * model_ntypes}, int_tensor_options).to(device);
        torch::Tensor dR_neigh_tensor = torch::from_blob(dR_neigh.data(), {1, inum, max_neighbor * model_ntypes, 4}, float_tensor_options).to(device,dtype);
        torch::Tensor atom_type_tensor = torch::from_blob(atom_types.data(), {ntypes}, int_tensor_options).to(device);
            
        torch::Tensor neighbor_type_list_tensor;
        if (model_name == "NEP") {
            // std::cout<< "=========do neighbor_type_list to tensor ==========="<<std::endl;
            neighbor_type_list_tensor = torch::from_blob(neighbor_type_list.data(), {1, inum, max_neighbor * model_ntypes}, int_tensor_options).to(device);
        }
        // auto t6 = std::chrono::high_resolution_clock::now();
        /*
        do forward for 4 models
        only 1 is used for MD
        1, 2, 3, 4 all for the deviation
        */

        for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
            if (ff_idx > 0 && (current_timestep % out_freq != 0)) continue;
            torch::Tensor Force;
            torch::Tensor Ei;
            torch::Tensor Etot;
            torch::Tensor Virial;
            if (model_name == "DP") {
                auto output = modules[ff_idx].forward({neighbor_list_tensor, imagetype_map_tensor, atom_type_tensor, dR_neigh_tensor, nghost}).toTuple();
                Force = output->elements()[2].toTensor().to(torch::kCPU);
                Ei = output->elements()[1].toTensor().to(torch::kCPU);
                if (ff_idx == 0) {
                    Etot = output->elements()[0].toTensor().to(torch::kCPU);
                    Virial = output->elements()[4].toTensor().to(torch::kCPU);            
                }
            } else if (model_name == "NEP") {
                // std::cout<< "=========do nep forward=========== nghost "<< nghost << std::endl;
                auto output = modules[ff_idx].forward({neighbor_list_tensor, imagetype_map_tensor, atom_type_tensor, dR_neigh_tensor, neighbor_type_list_tensor, nghost}).toTuple();
                Force = output->elements()[2].toTensor().to(torch::kCPU);
                Ei = output->elements()[1].toTensor().to(torch::kCPU);
                if (ff_idx == 0) {
                    Etot = output->elements()[0].toTensor().to(torch::kCPU);
                    Virial = output->elements()[4].toTensor().to(torch::kCPU);            
                }
            }
            auto F_ptr = Force.accessor<double, 3>();
            auto Ei_ptr = Ei.accessor<double, 2>();

            if (num_ff > 1 && (current_timestep % out_freq == 0)) {
                for (int i = 0; i < inum + nghost; i++)
                {
                    // f_n[ff_idx][i][0] = Force[0][i][0].item<double>();
                    // f_n[ff_idx][i][1] = Force[0][i][1].item<double>();
                    // f_n[ff_idx][i][2] = Force[0][i][2].item<double>();
                    f_n[ff_idx][i][0] = F_ptr[0][i][0];
                    f_n[ff_idx][i][1] = F_ptr[0][i][1];
                    f_n[ff_idx][i][2] = F_ptr[0][i][2];
                }
                for (int ii = 0; ii < inum; ii++) {
                    // e_atom_n[ff_idx][ii] = Ei[0][ii].item<double>();
                    e_atom_n[ff_idx][ii] = Ei_ptr[0][ii];
                }
            }

            if (ff_idx == 0) {
                // Etot = output->elements()[0].toTensor().to(torch::kCPU);
                // Virial = output->elements()[4].toTensor().to(torch::kCPU);
                // if (output->elements()[4].isTensor()) {
                //     calc_virial_from_mlff = true;
                //     torch::Tensor Virial = output->elements()[4].toTensor().to(torch::kCPU);
                // } else
                //     auto Virial = output->elements()[4];
                // get force

                // auto F_ptr = Force.accessor<double, 3>();
                // auto Ei_ptr = Ei.accessor<double, 2>();
                auto V_ptr = Virial.accessor<double, 2>();

                for (int i = 0; i < inum + nghost; i++)
                {
                    f[i][0] = F_ptr[0][i][0];
                    f[i][1] = F_ptr[0][i][1];
                    f[i][2] = F_ptr[0][i][2];
                }

                virial[0] = V_ptr[0][0];    // xx
                virial[1] = V_ptr[0][4];    // yy
                virial[2] = V_ptr[0][8];    // zz
                virial[3] = V_ptr[0][1];    // xy
                virial[4] = V_ptr[0][2];    // xz
                virial[5] = V_ptr[0][5];    // yz

                // get energy
                if (eflag) eng_vdwl = Etot[0][0].item<double>();

                if (eflag_atom)
                {
                    for (int ii = 0; ii < inum; ii++) {
                        eatom[ii] = Ei_ptr[0][ii];
                    }
                }
                // If virial needed calculate via F dot r.
                // if (vflag_fdotr) virial_fdotr_compute();
            }
        }
    } // if model_type == 0
    else if (model_type == 1) {
        double total_potential = 0.0;
        double total_virial[6] = {0.0};
        double* per_atom_potential = nullptr;
        double** per_atom_virial = nullptr;
        if (cvflag_atom) {
            per_atom_virial = cvatom;
        }
        if (eflag_atom) {
            per_atom_potential = eatom;
        }
        for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
            if ((num_ff == 1) or (current_timestep % out_freq != 0)) {
                                // can not set the atom->type (the type set in config) to nep forcefild order, because the ghost atoms type same as the conifg
                // The atomic types corresponding to the index of neighbors are constantly changing
                nep_cpu_models[ff_idx].compute_for_lammps(
                atom->nlocal, list->inum, list->ilist, list->numneigh, list->firstneigh, atom->type, atom->x,
                total_potential, total_virial, per_atom_potential, atom->f, per_atom_virial, ff_idx);
                                if (eflag) {
                    eng_vdwl += total_potential;
                }
                if (vflag) {
                    for (int component = 0; component < 6; ++component) {
                    virial[component] += total_virial[component];
                    }
                }
            } else {
                total_potential = 0.0;
                total_virial[6] = {0.0};
                // for multi models, the output step, should calculate deviation
                // 这里可以优化，快速置零
                for (int i = 0; i < list->inum + nghost; i++) {
                    f_n[ff_idx][i][0] = 0;
                    f_n[ff_idx][i][1] = 0;
                    f_n[ff_idx][i][2] = 0;
                    e_atom_n[ff_idx][i] = 0;
                }
                nep_cpu_models[ff_idx].compute_for_lammps(
                    atom->nlocal, list->inum, list->ilist, list->numneigh, list->firstneigh,atom->type, atom->x,
                    total_potential, total_virial, e_atom_n[ff_idx], f_n[ff_idx], per_atom_virial, ff_idx);
                if (ff_idx == 0) {
                    for (int i = 0; i < list->inum + nghost; i++) {
                        atom->f[i][0] = f_n[0][i][0];
                        atom->f[i][1] = f_n[0][i][1];
                        atom->f[i][2] = f_n[0][i][2];
                    }
                    if (eflag_atom) {
                        for (int i = 0; i < list->inum + nghost; i++) {
                           per_atom_potential[i] = e_atom_n[0][i];
                        }
                    }
                    if (eflag) {
                        eng_vdwl += total_potential;
                    }
                    if (vflag) {
                        for (int component = 0; component < 6; ++component) {
                            virial[component] += total_virial[component];
                        }
                    }    
                }
            } // else multi models out steps
        }   // for ff_idx      
    } // model_type == 1: nep_cpu version
    //   exploration mode.
    //   calculate the error of the force
    else if (model_type == 2) {
        double total_potential = 0.0;
        // double total_virial[6] = {0.0};
        double* per_atom_potential = nullptr;
        double** per_atom_virial = nullptr;
        if (cvflag_atom) {
            per_atom_virial = cvatom;
        }
        if (eflag_atom) {
            per_atom_potential = eatom;
        }
        
        // can not set the atom->type (the type set in config) to nep forcefild order, because the ghost atoms type same as the conifg
        // The atomic types corresponding to the index of neighbors are constantly changing

        std::vector<double> cpu_potential_per_atom(list->inum, 0.0);
        std::vector<double> cpu_force_per_atom(n_all * 3, 0.0);
        std::vector<double> cpu_total_virial(6, 0.0);

        // std::tie(itype_convert_map, neighbor_list, neighbor_angular_list, neigbor_num_list, neigbor_angular_num_list, rij_nep_gpu) = generate_neighdata_nep_gpu();
        // // for(int tmpi=0;tmpi< 10;tmpi++) {
        // //     printf("before ei [%d] = %f", tmpi, cpu_potential_per_atom[tmpi]);
        // // }
        // printf("before compute nall %d nlocal %d nghost %d inum %d\n", n_all, atom->nlocal, atom->nghost, list->inum);

        // nep_gpu_model.compute_small_box(
        // n_all, atom->nlocal, list->inum, nep_gpu_nm, itype_convert_map.data(), neigbor_num_list.data(), neighbor_list.data(), neigbor_angular_num_list.data(), neighbor_angular_list.data(), rij_nep_gpu.data(), 
        // cpu_potential_per_atom.data(), cpu_force_per_atom.data(), cpu_total_virial.data());

        std::tie(itype_convert_map, firstneighbor_cpu, position_cpu) = convert_dim(is_build_neighbor);
        // nep_gpu_model.compute_small_box_optim(
        // is_build_neighbor,
        // n_all, 
        // atom->nlocal,
        // list->inum, 
        // nep_gpu_nm, 
        // itype_convert_map.data(),
        // list->ilist,
        // list->numneigh,
        // firstneighbor_cpu.data(),
        // position_cpu.data(),
        // cpu_potential_per_atom.data(), 
        // cpu_force_per_atom.data(), 
        // cpu_total_virial.data());

        nep_gpu_model.compute_large_box_optim(
        is_build_neighbor,
        n_all, 
        atom->nlocal,
        list->inum, 
        nep_gpu_nm, 
        itype_convert_map.data(),
        list->ilist,
        list->numneigh,
        firstneighbor_cpu.data(),
        position_cpu.data(),
        cpu_potential_per_atom.data(), 
        cpu_force_per_atom.data(), 
        cpu_total_virial.data());


        // for(int tmpi=0;tmpi< 10;tmpi++) {
        //     printf("after ei [%d] = %f", tmpi, cpu_potential_per_atom[tmpi]);
        // }
        if (eflag) {
            // printf("before eng %f\n", eng_vdwl);
            double tmp_etot = 0;
            for (int i = 0; i < list->inum; ++i) {
            eng_vdwl += cpu_potential_per_atom[i];
            tmp_etot += cpu_potential_per_atom[i];
            }
            // printf("after eng %f real etot %f\n", eng_vdwl, tmp_etot);
            // std::cout<< "etot " << eng_vdwl << std::endl;
        }
        if (vflag) {
            for (int component = 0; component < 6; ++component) {
            virial[component] += cpu_total_virial[component];
            }
        }
        if (eflag_atom) {
            for (int i = 0; i < list->inum; ++i) {
                per_atom_potential[i] = cpu_potential_per_atom[i];
            }
        }
        // copy force
        for (int i = 0; i < n_all; ++i) {
            atom->f[i][0] = cpu_force_per_atom[i];
            atom->f[i][1] = cpu_force_per_atom[n_all + i];
            atom->f[i][2] = cpu_force_per_atom[2*n_all + i];

            // if (i % 1000 == 0) {
            //     printf("force_%d = [%f, %f, %f]\n", i, atom->f[i][0], atom->f[i][1], atom->f[i][2]);
            // }
        }
    }

    // for jit model or nep_cpu model, the mulit nep_cpu models have error
    if (num_ff > 1 && (current_timestep % out_freq == 0)) {
        // calculate model deviation with Force
        std::pair<double, double> result = calc_max_error(f_n, e_atom_n);
        max_err = result.first;
        max_err_ei = result.second;
        MPI_Allreduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, world);
        MPI_Allreduce(&max_err_ei, &global_max_err_ei, 1, MPI_DOUBLE, MPI_MAX, world);

        max_err_list.push_back(global_max_err);
        max_err_ei_list.push_back(global_max_err_ei);

        if (current_timestep % out_freq == 0) {
            if (me == 0) {
                fprintf(explrError_fp, "%9d %16.9f %16.9f\n", (max_err_list.size()-1)*out_freq, global_max_err, global_max_err_ei);
                fflush(explrError_fp);
            } 
        }
    } 
    // std::cout << "t4 " << (t5 - t4).count() * 0.000001 << "\tms" << std::endl;
    // std::cout << "t5 " << (t6 - t5).count() * 0.000001 << "\tms" << std::endl;
    // std::cout << "t6 " << (t7 - t6).count() * 0.000001 << "\tms" << std::endl;
}

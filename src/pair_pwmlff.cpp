#include <torch/script.h>
#include <torch/torch.h>
#include <iostream>
#include <chrono>
#include <cmath>
#include <cstring>

#include "pair_pwmlff.h"

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


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPWMLFF::PairPWMLFF(LAMMPS *lmp) : Pair(lmp)
{
    me = comm->me;
	writedata = 1;

}

PairPWMLFF::~PairPWMLFF()
{
    if (allocated) 
    {
        memory->destroy(setflag);
        memory->destroy(cutsq);
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
    if ((device_type == torch::kCUDA) && (rank == 0)) {
        if (num_devices < comm->nprocs) {
        std::cout << "----------------------------------------------------------------------------------" << std::endl;
        std::cout << " Warning: There are " << num_devices << " GPUs available " << std::endl;
        std::cout << " But have " << comm->nprocs << " MPI processes, may result in poor performance!!!" << std::endl;
        std::cout << "----------------------------------------------------------------------------------" << std::endl;
        }
    }
    dtype = torch::kFloat64;
    if (me == 0) utils::logmesg(this -> lmp, "<---- Loading model ---->");
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        std::string model_file = models[ff_idx];
        try
        {
            // module = torch::jit::load(model_file, c10::Device(device));
            module = torch::jit::load(model_file, device);
            module.to(dtype);
            // module.eval();
            modules.push_back(module);
            if (me == 0) printf("\nLoading model file:   %s\n", model_file.c_str());
        }
        catch (const c10::Error e)
        {
            std::cerr << "Failed to load model :" << e.msg() << std::endl;
        }
    }
    cutoff = module.attr("Rmax").toDouble();
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
    if (force->newton_pair == 0) error->all(FLERR, "Pair style PWMATMLFF requires newon pair on");
    // Using a nearest neighbor table of type full
    neighbor->add_request(this, NeighConst::REQ_FULL);
}
/* ---------------------------------------------------------------------- */

std::pair<double, double> PairPWMLFF::calc_max_error(std::vector<torch::Tensor> all_forces, std::vector<torch::Tensor> all_ei)
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

    forces_accessors.clear();
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        forces_accessors.push_back(all_forces[ff_idx].accessor<double, 3>());
        // p_ff_idx is for reverse comm
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
        // auto F_ptr = all_forces[ff_idx].accessor<double, 3>();
        auto Ei_ptr = all_ei[ff_idx].accessor<double, 2>();
        for (i = 0; i < nlocal; i++) {
            // std::cout << "!!!forces_accessors[ff_idx][0][i][0] = " << tag[i] << " " << forces_accessors[ff_idx][0][i][0] << std::endl;
            f_ave[i * 3 + 0] += forces_accessors[ff_idx][0][i][0];
            f_ave[i * 3 + 1] += forces_accessors[ff_idx][0][i][1];
            f_ave[i * 3 + 2] += forces_accessors[ff_idx][0][i][2];
            ei_ave[i] += Ei_ptr[0][i];
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
        // auto F_ptr = all_forces[ff_idx].accessor<double, 3>();
        auto Ei_ptr = all_ei[ff_idx].accessor<double, 2>();
        for (i = 0; i < nlocal; i++) {
            // std::cout << "???forces_accessors[ff_idx][0][i][0] = " << tag[i] << " " << forces_accessors[ff_idx][0][i][0] << std::endl;
            f_err[ff_idx][i * 3 + 0] = forces_accessors[ff_idx][0][i][0] - f_ave[i * 3 + 0];
            f_err[ff_idx][i * 3 + 1] = forces_accessors[ff_idx][0][i][1] - f_ave[i * 3 + 1];
            f_err[ff_idx][i * 3 + 2] = forces_accessors[ff_idx][0][i][2] - f_ave[i * 3 + 2];
            ei_err[ff_idx][i] = Ei_ptr[0][i] - ei_ave[i];
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
    // auto F_ptr = all_forces[p_ff_idx].accessor<double, 3>();
    for (i = first; i < last; i++) {
        buf[m++] = forces_accessors[p_ff_idx][0][i][0];
        buf[m++] = forces_accessors[p_ff_idx][0][i][1];
        buf[m++] = forces_accessors[p_ff_idx][0][i][2];
    }
    return m;
}

void PairPWMLFF::unpack_reverse_comm(int n, int* list, double* buf) {
    int i, j, m;

    m = 0;
    // auto F_ptr = all_forces[p_ff_idx].accessor<double, 3>();
    for (i = 0; i < n; i++) {
        j = list[i];
        forces_accessors[p_ff_idx][0][j][0] += buf[m++];
        forces_accessors[p_ff_idx][0][j][1] += buf[m++];
        forces_accessors[p_ff_idx][0][j][2] += buf[m++];
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

void PairPWMLFF::compute(int eflag, int vflag)
{
    if (eflag || vflag) ev_setup(eflag, vflag);
    double **x = atom->x;
    double **f = atom->f;
    int *type = atom->type;
    int newton_pair = force->newton_pair;
    double *virial = force->pair->virial;
    int ff_idx;
    // int nlocal = atom->nlocal;
    int current_timestep = update->ntimestep;
    // int total_timestep = update->laststep;
    int ntypes = atom->ntypes;
    int nghost = atom->nghost;
    // int n_all = nlocal + nghost;
    bool calc_virial_from_mlff = false;
    bool calc_egroup_from_mlff = false;

    int inum, jnum, itype, jtype;
    double max_err, global_max_err, max_err_ei, global_max_err_ei;    
    inum = list->inum;

    // auto t4 = std::chrono::high_resolution_clock::now();
    auto [imagetype_map, neighborlist, dR_neigh] = generate_neighdata();
    if (inum == 0) return;
    // auto t5 = std::chrono::high_resolution_clock::now();
    auto int_tensor_options = torch::TensorOptions().dtype(torch::kInt);
    auto float_tensor_options = torch::TensorOptions().dtype(torch::kFloat64);
    torch::Tensor imagetype_map_tensor = torch::from_blob(imagetype_map.data(), {inum}, int_tensor_options).to(device);
    torch::Tensor neighbor_list_tensor = torch::from_blob(neighborlist.data(), {1, inum, max_neighbor * model_ntypes}, int_tensor_options).to(device);
    torch::Tensor dR_neigh_tensor = torch::from_blob(dR_neigh.data(), {1, inum, max_neighbor * model_ntypes, 4}, float_tensor_options).to(device,dtype);
    torch::Tensor atom_type_tensor = torch::from_blob(atom_types.data(), {ntypes}, int_tensor_options).to(device);
    // auto t6 = std::chrono::high_resolution_clock::now();
    /*
      do forward for 4 models
      1 is used for MD
      2, 3, 4 for the test
    */
    all_forces.clear();     // clear the force vector
    all_ei.clear();         // clear the atomic energy vector
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        auto output = modules[ff_idx].forward({neighbor_list_tensor, imagetype_map_tensor, atom_type_tensor, dR_neigh_tensor, nghost}).toTuple();
            torch::Tensor Force = output->elements()[2].toTensor().to(torch::kCPU);
            torch::Tensor Ei = output->elements()[1].toTensor().to(torch::kCPU);
            all_forces.push_back(Force);
            all_ei.push_back(Ei);
        if (ff_idx == 0) {
            torch::Tensor Etot = output->elements()[0].toTensor().to(torch::kCPU);
            torch::Tensor Virial = output->elements()[4].toTensor().to(torch::kCPU);
            // if (output->elements()[4].isTensor()) {
            //     calc_virial_from_mlff = true;
            //     torch::Tensor Virial = output->elements()[4].toTensor().to(torch::kCPU);
            // } else
            //     auto Virial = output->elements()[4];
            // get force
            auto F_ptr = Force.accessor<double, 3>();
            auto Ei_ptr = Ei.accessor<double, 2>();
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
    // auto t7 = std::chrono::high_resolution_clock::now();
    /*
      exploration mode.
      calculate the error of the force
  */
    if (num_ff > 1) {
        // calculate model deviation with Force
        std::pair<double, double> result = calc_max_error(all_forces, all_ei);
        max_err = result.first;
        max_err_ei = result.second;
        MPI_Allreduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, world);
        MPI_Allreduce(&max_err_ei, &global_max_err_ei, 1, MPI_DOUBLE, MPI_MAX, world);

        max_err_list.push_back(global_max_err);
        max_err_ei_list.push_back(global_max_err_ei);

        if (current_timestep % out_freq == 0) {
            if (me == 0) {
                fprintf(explrError_fp, "%9d %16.9f %16.9f\n", max_err_list.size()-1, global_max_err, global_max_err_ei);
                fflush(explrError_fp);
            } 
        }
    } 
    // std::cout << "t4 " << (t5 - t4).count() * 0.000001 << "\tms" << std::endl;
    // std::cout << "t5 " << (t6 - t5).count() * 0.000001 << "\tms" << std::endl;
    // std::cout << "t6 " << (t7 - t6).count() * 0.000001 << "\tms" << std::endl;

}

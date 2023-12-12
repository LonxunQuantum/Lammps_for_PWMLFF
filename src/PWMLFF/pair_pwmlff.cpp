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
    int iarg;

    if (narg <= 0) error->all(FLERR, "Illegal pair_style command");
    std::vector<std::string> models;

    while (iarg < narg) {
        iarg++;
    }
    num_ff = utils::inumeric(FLERR, arg[0], false, lmp);    // number of models
    for (int ii = 1; ii < iarg; ++ii) {
        models.push_back(arg[ii]);                          // model files
    }
    
    device = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;
    dtype = torch::kFloat64;
    if (me == 0) utils::logmesg(this -> lmp, "<---- Loading model ---->");
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        std::string model_file = models[ff_idx];
        try
        {
            module = torch::jit::load(model_file, c10::Device(device));
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
    utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

    int count = 0;
    for(int i = ilo; i <= ihi; i++) {
        for(int j = MAX(jlo,i); j <= jhi; j++) 
        {
            setflag[i][j] = 1;
            count++;
        }
    }

    auto type_map_module = module.attr("atom_type").toList();
    if (ntypes > type_map_module.size() || ntypes != narg - 2)
    {
        error->all(FLERR, "Element mapping is not correct, ntypes = " + std::to_string(ntypes));
    }
    for (int ii = 2; ii < narg; ++ii) {
        int temp = std::stoi(arg[ii]);
        auto iter = std::find(type_map_module.begin(), type_map_module.end(), temp);   
        if (iter != type_map_module.end() || arg[ii] == 0)
        {
            type_map.push_back(temp);
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
    // Using a nearest neighbor table of type full
    neighbor->add_request(this, NeighConst::REQ_FULL);
}
/* ---------------------------------------------------------------------- */

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
    int ntypes = atom->ntypes;
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
    use_type.resize(n_all);

    
    for (int ii = 0; ii < n_all; ii++)
    {
        use_type[ii] = type_map[type[ii] - 1];
        // type[0], type[1], type[2], type[3], type[4], : 2, 2, 1, 2, 2, ...
        // type_map[0], type_map[1] : 6, 1
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
        itype = type[i];                        // 2, 2, 1, 2, 2
        jlist = firstneigh[i];
        jnum = numneigh[i];
        imagetype_map[ii] = itype - 1;          // 1, 1, 0, 1, 1        python index from 0
        // imagetype[ii] = use_type[i];            // 1, 1, 6, 1, 1

        for (jj = 0; jj < jnum; jj++)
        {
            j = jlist[jj];                      // 1, 2, 3, 4;   0, 2, 3, 4;   0, 1, 3, 4;   0, 1, 2, 4;   0, 1, 2, 3
            delx = x[j][0] - x[i][0];
            dely = x[j][1] - x[i][1];
            delz = x[j][2] - x[i][2];
            rsq = delx * delx + dely * dely + delz * delz;
            jtype = type[j];                    // 
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
                if (rsq < min_dR) min_dR = rsq;
            }
        }
    }

    MPI_Allreduce(&min_dR, &min_dR_all, 1, MPI_DOUBLE, MPI_MIN, world);

    if (min_dR_all < 0.81) {
        if (me == 0) {
            std::cout << "ERROR: there are two atoms too close, min_dR_all = " << min_dR_all << std::endl;
            error->universe_all(FLERR, "there are two atoms too close");
        }
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
    int ntypes = atom->ntypes;
    int nghost = atom->nghost;
    // int n_all = nlocal + nghost;
    bool calc_virial_from_mlff = false;
    bool calc_egroup_from_mlff = false;

    int inum, jnum, itype, jtype;
        
    inum = list->inum;

    // auto t4 = std::chrono::high_resolution_clock::now();
    auto [imagetype_map, neighborlist, dR_neigh] = generate_neighdata();
    // auto t5 = std::chrono::high_resolution_clock::now();
    auto int_tensor_options = torch::TensorOptions().dtype(torch::kInt);
    auto float_tensor_options = torch::TensorOptions().dtype(torch::kFloat64);
    torch::Tensor imagetype_map_tensor = torch::from_blob(imagetype_map.data(), {inum}, int_tensor_options).to(device);
    torch::Tensor neighbor_list_tensor = torch::from_blob(neighborlist.data(), {1, inum, max_neighbor * ntypes}, int_tensor_options).to(device);
    torch::Tensor dR_neigh_tensor = torch::from_blob(dR_neigh.data(), {1, inum, max_neighbor * ntypes, 4}, float_tensor_options).to(device,dtype);
    torch::Tensor type_map_tensor = torch::from_blob(type_map.data(), {ntypes}, int_tensor_options).to(device);
    // auto t6 = std::chrono::high_resolution_clock::now();
    /*
      do forward for 4 models
      1 is used for MD
      2, 3, 4 for the test
    */
    for (ff_idx = 0; ff_idx < num_ff; ff_idx++) {
        auto output = modules[ff_idx].forward({neighbor_list_tensor, imagetype_map_tensor, type_map_tensor, dR_neigh_tensor, nghost}).toTuple();
        if (ff_idx == 0) {
            torch::Tensor Etot = output->elements()[0].toTensor().to(torch::kCPU);
            torch::Tensor Ei = output->elements()[1].toTensor().to(torch::kCPU);
            torch::Tensor Force = output->elements()[2].toTensor().to(torch::kCPU);
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
      select candidates
      write as fractional coordinate
      Note: only write data at the very end (in the destructor)!
  */
    // std::cout << "t4 " << (t5 - t4).count() * 0.000001 << "\tms" << std::endl;
    // std::cout << "t5 " << (t6 - t5).count() * 0.000001 << "\tms" << std::endl;
    // std::cout << "t6 " << (t7 - t6).count() * 0.000001 << "\tms" << std::endl;

}

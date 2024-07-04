/*
    Copyright 2017 Zheyong Fan, Ville Vierimaa, Mikko Ervasti, and Ari Harju
    This file is part of GPUMD.
    GPUMD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    GPUMD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with GPUMD.  If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------80
The neuroevolution potential (NEP)
Ref: Zheyong Fan et al., Neuroevolution machine learning potentials:
Combining high accuracy and low cost in atomistic simulations and application to
heat transport, Phys. Rev. B. 104, 104309 (2021).
------------------------------------------------------------------------------*/

#include "nep3.cuh"
#include "nep3_small_box.cuh"
#include "utilities/common.cuh"
#include "utilities/error.cuh"
#include "utilities/nep_utilities.cuh"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

const std::string ELEMENTS[NUM_ELEMENTS] = {
  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",
  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
  "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re",
  "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
  "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"};

void checkArray(const char* name, void* ptr, size_t size) {
    printf("%s: pointer = %p, size = %zu\n", name, ptr, size);
}

NEP3::NEP3() {}

void NEP3::init_from_file(const char* file_potential, const bool is_rank_0, const int in_device_id)
{
  rank_0 = is_rank_0;
  device_id = in_device_id;
  atom_nums = 0;
  std::ifstream input(file_potential);
  if (!input.is_open()) {
    std::cout << "Failed to open " << file_potential << std::endl;
    exit(1);
  }

  // nep3 1 C
  std::vector<std::string> tokens = get_tokens(input);
  if (tokens.size() < 3) {
    std::cout << "The first line of nep.txt should have at least 3 items." << std::endl;
    exit(1);
  }
  if (tokens[0] == "nep4") {
    paramb.version = 4;
  }
  paramb.num_types = get_int_from_token(tokens[1], __FILE__, __LINE__);
  if (tokens.size() != 2 + paramb.num_types) {
    std::cout << "The first line of nep.txt should have " << paramb.num_types << " atom symbols."
              << std::endl;
    exit(1);
  }

  if (paramb.num_types == 1) {
    printf("Use the NEP%d potential with %d atom type.\n", paramb.version, paramb.num_types);
  } else {
    printf("Use the NEP%d potential with %d atom types.\n", paramb.version, paramb.num_types);
  }
  element_atomic_number_list.resize(paramb.num_types);
  for (int n = 0; n < paramb.num_types; ++n) {
    int atomic_number = 0;
    for (int m = 0; m < NUM_ELEMENTS; ++m) {
      if (tokens[2 + n] == ELEMENTS[m]) {
        atomic_number = m + 1;
        break;
      }
    }
    element_atomic_number_list[n] = atomic_number;
    printf("    type %d (%s).\n", n, tokens[2 + n].c_str());
  }

  // cutoff 4.2 3.7 80 47
  tokens = get_tokens(input);
  if (tokens.size() != 3 && tokens.size() != 5) {
    std::cout << "This line should be cutoff rc_radial rc_angular [MN_radial] [MN_angular].\n";
    exit(1);
  }
  paramb.rc_radial = get_float_from_token(tokens[1], __FILE__, __LINE__);
  paramb.rc_angular = get_float_from_token(tokens[2], __FILE__, __LINE__);
  printf("    radial cutoff = %g A.\n", paramb.rc_radial);
  printf("    angular cutoff = %g A.\n", paramb.rc_angular);

  paramb.MN_radial = 500;
  paramb.MN_angular = 100;

  if (tokens.size() == 5) {
    int MN_radial = get_int_from_token(tokens[3], __FILE__, __LINE__);
    int MN_angular = get_int_from_token(tokens[4], __FILE__, __LINE__);
    printf("    MN_radial = %d.\n", MN_radial);
    printf("    MN_angular = %d.\n", MN_angular);
    paramb.MN_radial = int(ceil(MN_radial * 1.25));
    paramb.MN_angular = int(ceil(MN_angular * 1.25));
    printf("    enlarged MN_radial = %d.\n", paramb.MN_radial);
    printf("    enlarged MN_angular = %d.\n", paramb.MN_angular);
  }

  // n_max 10 8
  tokens = get_tokens(input);
  if (tokens.size() != 3) {
    std::cout << "This line should be n_max n_max_radial n_max_angular." << std::endl;
    exit(1);
  }
  paramb.n_max_radial = get_int_from_token(tokens[1], __FILE__, __LINE__);
  paramb.n_max_angular = get_int_from_token(tokens[2], __FILE__, __LINE__);
  printf("    n_max_radial = %d.\n", paramb.n_max_radial);
  printf("    n_max_angular = %d.\n", paramb.n_max_angular);

  // basis_size 10 8
  if (paramb.version >= 3) {
    tokens = get_tokens(input);
    if (tokens.size() != 3) {
      std::cout << "This line should be basis_size basis_size_radial basis_size_angular."
                << std::endl;
      exit(1);
    }
    paramb.basis_size_radial = get_int_from_token(tokens[1], __FILE__, __LINE__);
    paramb.basis_size_angular = get_int_from_token(tokens[2], __FILE__, __LINE__);
    printf("    basis_size_radial = %d.\n", paramb.basis_size_radial);
    printf("    basis_size_angular = %d.\n", paramb.basis_size_angular);
  }

  // l_max
  tokens = get_tokens(input);
  if (paramb.version == 2) {
    if (tokens.size() != 2) {
      std::cout << "This line should be l_max l_max_3body." << std::endl;
      exit(1);
    }
  } else {
    if (tokens.size() != 4) {
      std::cout << "This line should be l_max l_max_3body l_max_4body l_max_5body." << std::endl;
      exit(1);
    }
  }

  paramb.L_max = get_int_from_token(tokens[1], __FILE__, __LINE__);
  printf("    l_max_3body = %d.\n", paramb.L_max);
  paramb.num_L = paramb.L_max;

  if (paramb.version >= 3) {
    int L_max_4body = get_int_from_token(tokens[2], __FILE__, __LINE__);
    int L_max_5body = get_int_from_token(tokens[3], __FILE__, __LINE__);
    printf("    l_max_4body = %d.\n", L_max_4body);
    printf("    l_max_5body = %d.\n", L_max_5body);
    if (L_max_4body == 2) {
      paramb.num_L += 1;
    }
    if (L_max_5body == 1) {
      paramb.num_L += 1;
    }
  }

  paramb.dim_angular = (paramb.n_max_angular + 1) * paramb.num_L;

  // ANN
  tokens = get_tokens(input);
  if (tokens.size() != 3) {
    std::cout << "This line should be ANN num_neurons 0." << std::endl;
    exit(1);
  }
  annmb.num_neurons1 = get_int_from_token(tokens[1], __FILE__, __LINE__);
  annmb.dim = (paramb.n_max_radial + 1) + paramb.dim_angular;
  if (paramb.model_type == 3) {
    annmb.dim += 1;
  }
  printf("    ANN = %d-%d-1.\n", annmb.dim, annmb.num_neurons1);

  // calculated parameters:
  // rc = paramb.rc_radial; // largest cutoff
  paramb.rcinv_radial = 1.0f / paramb.rc_radial;
  paramb.rcinv_angular = 1.0f / paramb.rc_angular;
  paramb.num_types_sq = paramb.num_types * paramb.num_types;
  annmb.num_para =
    (annmb.dim + 2) * annmb.num_neurons1 * (paramb.version == 4 ? paramb.num_types : 1) + (paramb.version == 4 ? paramb.num_types : 1);
  if (paramb.model_type == 2) {
    // Polarizability models have twice as many parameters
    annmb.num_para *= 2;
  }
  printf("    number of neural network parameters = %d.\n", annmb.num_para);
  annmb.num_c2   = paramb.num_types_sq * (paramb.n_max_radial + 1) * (paramb.basis_size_radial + 1);
  annmb.num_c3   = paramb.num_types_sq * (paramb.n_max_angular + 1) * (paramb.basis_size_angular + 1);
  int num_para_descriptor =annmb.num_c2 + annmb.num_c3;
    // paramb.num_types_sq * ((paramb.n_max_radial + 1) * (paramb.basis_size_radial + 1) +
    //                        (paramb.n_max_angular + 1) * (paramb.basis_size_angular + 1));

  printf("    number of descriptor parameters = %d.\n", num_para_descriptor);
  annmb.num_para += num_para_descriptor;
  printf("    total number of parameters = %d.\n", annmb.num_para);

  paramb.num_c_radial =
    paramb.num_types_sq * (paramb.n_max_radial + 1) * (paramb.basis_size_radial + 1);

  // NN and descriptor parameters
  std::vector<float> parameters(annmb.num_para);
  for (int n = 0; n < annmb.num_para; ++n) {
    tokens = get_tokens(input);
    parameters[n] = get_float_from_token(tokens[0], __FILE__, __LINE__);
    // std::cout<<"nn param " << n << " " << parameters[n] << std::endl;
  }
  nep_data.parameters.resize(annmb.num_para);
  nep_data.parameters.copy_from_host(parameters.data());
  update_potential(nep_data.parameters.data(), annmb);

  for (int d = 0; d < annmb.dim; ++d) {
    tokens = get_tokens(input);
    paramb.q_scaler[d] = get_float_from_token(tokens[0], __FILE__, __LINE__);
    // std::cout<<"q_scaler " << d << " " << paramb.q_scaler[d] << std::endl;
  }

  // nep_data.f12x.resize(num_atoms * paramb.MN_angular);
  // nep_data.f12y.resize(num_atoms * paramb.MN_angular);
  // nep_data.f12z.resize(num_atoms * paramb.MN_angular);
  
  // nep_data.NN_radial.resize(num_atoms);
  // nep_data.NL_radial.resize(num_atoms * paramb.MN_radial);
  // nep_data.NN_angular.resize(num_atoms);
  // nep_data.NL_angular.resize(num_atoms * paramb.MN_angular);
  // nep_data.Fp.resize(num_atoms * annmb.dim);
  // nep_data.sum_fxyz.resize(num_atoms * (paramb.n_max_angular + 1) * NUM_OF_ABC);
  
  // nep_data.cell_count.resize(num_atoms);
  // nep_data.cell_count_sum.resize(num_atoms);
  // nep_data.cell_contents.resize(num_atoms);
  // nep_data.cpu_NN_radial.resize(num_atoms);
  // nep_data.cpu_NN_angular.resize(num_atoms);
#ifdef USE_TABLE
  construct_table(parameters.data());
  printf("    use tabulated radial functions to speed up.\n");
#endif
}

NEP3::~NEP3(void)
{
  // nothing
}

// void NEP3:: set_partition(int n_all, int n_local, int n_ghost, int n_inum) {
//   nlocal = n_local;
//   nghost = n_ghost;
//   inum = n_inum;
//   if (inum > nlocal) {
//     partition = inum;
//   }
//   else {
//     partition = nlocal;
//   }
//   nall = n_all;
// }

void NEP3::checkMemoryUsage(int sgin) {
  // if (rank_0) {
    size_t free_mem, total_mem;
    cudaError_t error = cudaMemGetInfo(&free_mem, &total_mem);
    if (error != cudaSuccess) {
        std::cerr << "cudaMemGetInfo failed: " << cudaGetErrorString(error) << std::endl;
        return;
    }
    std::cout << device_id << " Free memory: "  << sgin << " " << free_mem / (1024.0 * 1024.0) << " MB" << std::endl;
    // std::cout << device_id << " Total memory: " << sgin << " " << total_mem / (1024.0 * 1024.0) << " MB" << std::endl;
  // }
}

void NEP3::rest_nep_data(int input_atom_num, int n_local, int n_all, int max_neighbor, bool build_neighbor) {
  // printf("(build neighbor %d) (atom_nums %d change %d) (atom_local %d change %d) (atom_nums_all %d change %d)\n", build_neighbor, input_atom_num, (atom_nums != input_atom_num), n_local, (atom_nlocal != n_local), n_all, (atom_nums_all != n_all));
  if (atom_nums != input_atom_num) {
    atom_nums = input_atom_num;
    nep_data.NN_radial.resize(atom_nums);
    nep_data.NL_radial.resize(atom_nums * paramb.MN_radial);
    nep_data.NN_angular.resize(atom_nums);
    nep_data.NL_angular.resize(atom_nums * paramb.MN_angular);
    nep_data.potential_per_atom.resize(atom_nums);
    lmp_data.ilist.resize(atom_nums);
    lmp_data.numneigh.resize(atom_nums);
    lmp_data.firstneigh.resize(atom_nums*max_neighbor);
    lmp_data.r12.resize(atom_nums*max_neighbor*6);
  }

  if (atom_nlocal != n_local) {
    atom_nlocal = n_local;
    nep_data.Fp.resize(n_local * annmb.dim);
    nep_data.sum_fxyz.resize(n_local * (paramb.n_max_angular + 1) * NUM_OF_ABC);
  }

  if (atom_nums_all != n_all) {
    atom_nums_all = n_all;
    nep_data.force_per_atom.resize(atom_nums_all * 3);
    nep_data.virial_per_atom.resize(atom_nums_all * 9);
    nep_data.total_virial.resize(6);
    lmp_data.type.resize(atom_nums_all);
    lmp_data.position.resize(atom_nums_all*3);
  }
  nep_data.potential_per_atom.fill(0.0);
  nep_data.force_per_atom.fill(0.0);
  nep_data.virial_per_atom.fill(0.0);
  nep_data.total_virial.fill(0.0);
}

void NEP3::update_potential(float* parameters, ANN& ann)
{
  float* pointer = parameters;
  for (int t = 0; t < paramb.num_types; ++t) {
    if (t > 0 && paramb.version != 4) { // Use the same set of NN parameters for NEP2 and NEP3
      pointer -= (ann.dim + 2) * ann.num_neurons1;
    }
    ann.w0[t] = pointer;
    pointer += ann.num_neurons1 * ann.dim;
    ann.b0[t] = pointer;
    pointer += ann.num_neurons1;
    ann.w1[t] = pointer;
    pointer += ann.num_neurons1;
  }
  ann.b1 = pointer;
  // pointer += 1;
  pointer += (paramb.version == 4 ? paramb.num_types : 1);

  ann.c = pointer;
}

#ifdef USE_TABLE
void NEP3::construct_table(float* parameters)
{
  nep_data.gn_radial.resize(table_length * paramb.num_types_sq * (paramb.n_max_radial + 1));
  nep_data.gnp_radial.resize(table_length * paramb.num_types_sq * (paramb.n_max_radial + 1));
  nep_data.gn_angular.resize(table_length * paramb.num_types_sq * (paramb.n_max_angular + 1));
  nep_data.gnp_angular.resize(table_length * paramb.num_types_sq * (paramb.n_max_angular + 1));
  std::vector<float> gn_radial(table_length * paramb.num_types_sq * (paramb.n_max_radial + 1));
  std::vector<float> gnp_radial(table_length * paramb.num_types_sq * (paramb.n_max_radial + 1));
  std::vector<float> gn_angular(table_length * paramb.num_types_sq * (paramb.n_max_angular + 1));
  std::vector<float> gnp_angular(table_length * paramb.num_types_sq * (paramb.n_max_angular + 1));
  float* c_pointer =
    parameters +
    (annmb.dim + 2) * annmb.num_neurons1 * (paramb.version == 4 ? paramb.num_types : 1) + 1;
  construct_table_radial_or_angular(
    paramb.version,
    paramb.num_types,
    paramb.num_types_sq,
    paramb.n_max_radial,
    paramb.basis_size_radial,
    paramb.rc_radial,
    paramb.rcinv_radial,
    c_pointer,
    gn_radial.data(),
    gnp_radial.data());
  construct_table_radial_or_angular(
    paramb.version,
    paramb.num_types,
    paramb.num_types_sq,
    paramb.n_max_angular,
    paramb.basis_size_angular,
    paramb.rc_angular,
    paramb.rcinv_angular,
    c_pointer + paramb.num_c_radial,
    gn_angular.data(),
    gnp_angular.data());
  nep_data.gn_radial.copy_from_host(gn_radial.data());
  nep_data.gnp_radial.copy_from_host(gnp_radial.data());
  nep_data.gn_angular.copy_from_host(gn_angular.data());
  nep_data.gnp_angular.copy_from_host(gnp_angular.data());
}
#endif

// small box possibly used for active learning:
void NEP3::compute_small_box(
  int n_all, //n_local + nghost
  int nlocal, 
  int N, // list->inum
  int NM,// maxneighbors
  int* itype,//atoms' type,the len is n_all
  int* cpu_nn_radail, // the len is N, value is the neighbor nums of atom_i 
  int* cpu_nl_radail, // the len is N*NM, value is the neighbor list of atom_i 
  int* cpu_nn_angular,
  int* cpu_nl_angular,
  float* cpu_r12,   // the len is N*NM*6, value is the neighbor list of atom_i 
  double* cpu_potential_per_atom, // the output of ei
  double* cpu_force_per_atom,     // the output of force
  double* cpu_total_virial     // the output of virial len 6
  )
{
//   rest_nep_data(N, nlocal, n_all, NM);

//   int N1 = 0;
//   const int BLOCK_SIZE = 64;
//   // const int N = type.size();
//   const int grid_size = (N - 1) / BLOCK_SIZE + 1;

//   // const int big_neighbor_size = 2000;
//   // const int size_x12 = type.size() * big_neighbor_size;
//   const int size_x12 = N * NM;
//   GPU_Vector<float> r12(size_x12 * 6);
//   GPU_Vector<int> type(n_all);

//   nep_data.NN_radial.copy_from_host(cpu_nn_radail);
//   nep_data.NL_radial.copy_from_host(cpu_nl_radail);
//   nep_data.NN_angular.copy_from_host(cpu_nn_angular);
//   nep_data.NL_angular.copy_from_host(cpu_nl_angular);
//   r12.copy_from_host(cpu_r12);
//   type.copy_from_host(itype);

//   GPU_Vector<double> potential_per_atom(N);
//   potential_per_atom.fill(0.0);
//   GPU_Vector<double> force_per_atom(n_all * 3);
//   force_per_atom.fill(0.0);
//   GPU_Vector<double> virial_per_atom(n_all * N * 9);
//   virial_per_atom.fill(0.0);
//   GPU_Vector<double> total_virial(6);
//   total_virial.fill(0.0);

//   // printf("before find_descriptor_small_box\n");
//   find_descriptor_small_box<<<grid_size, BLOCK_SIZE>>>(
//     paramb,
//     annmb,
//     N,
//     N1,
//     nlocal,
//     nep_data.NN_radial.data(),
//     nep_data.NL_radial.data(),
//     nep_data.NN_angular.data(),
//     nep_data.NL_angular.data(),
//     type.data(),
//     r12.data(),
//     r12.data() + size_x12,
//     r12.data() + size_x12 * 2,
//     r12.data() + size_x12 * 3,
//     r12.data() + size_x12 * 4,
//     r12.data() + size_x12 * 5,
//     false,//is_polarizability
// #ifdef USE_TABLE
//     nep_data.gn_radial.data(),
//     nep_data.gn_angular.data(),
// #endif
//     potential_per_atom.data(),
//     nep_data.Fp.data(),
//     virial_per_atom.data(),
//     nep_data.sum_fxyz.data());
//   CUDA_CHECK_KERNEL
//   // bool is_dipole = paramb.model_type == 1;
//   find_force_radial_small_box<<<grid_size, BLOCK_SIZE>>>(
//     paramb,
//     annmb,
//     N,
//     N1,
//     nlocal,
//     nep_data.NN_radial.data(),
//     nep_data.NL_radial.data(),
//     type.data(),
//     r12.data(),
//     r12.data() + size_x12,
//     r12.data() + size_x12 * 2,
//     nep_data.Fp.data(),
//     false, //is_dipole,
// #ifdef USE_TABLE
//     nep_data.gnp_radial.data(),
// #endif
//     force_per_atom.data(),
//     force_per_atom.data() + n_all,
//     force_per_atom.data() + n_all * 2,
//     virial_per_atom.data(),
//     total_virial.data());
//   CUDA_CHECK_KERNEL
//   find_force_angular_small_box<<<grid_size, BLOCK_SIZE>>>(
//     paramb,
//     annmb,
//     N,
//     N1,
//     nlocal,
//     nep_data.NN_angular.data(),
//     nep_data.NL_angular.data(),
//     type.data(),
//     r12.data() + size_x12 * 3,
//     r12.data() + size_x12 * 4,
//     r12.data() + size_x12 * 5,
//     nep_data.Fp.data(),
//     nep_data.sum_fxyz.data(),
//     false,//is_dipole,
// #ifdef USE_TABLE
//     nep_data.gn_angular.data(),
//     nep_data.gnp_angular.data(),
// #endif
//     force_per_atom.data(),
//     force_per_atom.data() + n_all,
//     force_per_atom.data() + n_all * 2,
//     virial_per_atom.data(),
//     total_virial.data());
//   CUDA_CHECK_KERNEL
  
//   total_virial.copy_to_host(cpu_total_virial);
//   potential_per_atom.copy_to_host(cpu_potential_per_atom);
//   force_per_atom.copy_to_host(cpu_force_per_atom);
}



// small box possibly used for active learning:
void NEP3::compute_small_box_optim(
  bool is_build_neighbor,
  int n_all, //n_local + nghost
  int nlocal,
  int N, //atom nums
  int NM,// maxneighbors
  int* itype_cpu,//atoms' type,the len is [n_all]
  int* ilist_cpu, // atom i list
  int* numneigh_cpu, // the neighbor nums of each i, [inum]
  int* firstneigh_cpu, // the neighbor list of each i, [inum * NM]
  double* position_cpu, // postion of atoms x, [n_all * 3]
  double* cpu_potential_per_atom, // the output of ei
  double* cpu_force_per_atom,     // the output of force
  double* cpu_total_virial     // the output of virial
) {
  int N1 = 0;
  // checkMemoryUsage(0);
  rest_nep_data(N, nlocal, n_all, NM, is_build_neighbor);

  const int BLOCK_SIZE = 64;
  const int grid_size = (N - 1) / BLOCK_SIZE + 1;
  const int size_x12 = N * NM;

  if (is_build_neighbor) {
    lmp_data.type.copy_from_host(itype_cpu);
    lmp_data.ilist.copy_from_host(ilist_cpu);
    lmp_data.numneigh.copy_from_host(numneigh_cpu);
    lmp_data.firstneigh.copy_from_host(firstneigh_cpu);
    // checkMemoryUsage(2);
  }
  lmp_data.position.copy_from_host(position_cpu);
  
  find_neighbor_list_small_box<<<grid_size, BLOCK_SIZE>>>(
    paramb,
    n_all,
    N,
    N1,
    NM,
    lmp_data.type.data(),
    lmp_data.ilist.data(),
    lmp_data.numneigh.data(),
    lmp_data.firstneigh.data(),
    lmp_data.position.data(),
    lmp_data.position.data() + n_all,
    lmp_data.position.data() + n_all * 2,
    nep_data.NN_radial.data(),
    nep_data.NL_radial.data(),
    nep_data.NN_angular.data(),
    nep_data.NL_angular.data(),
    lmp_data.r12.data(),
    lmp_data.r12.data() + size_x12,
    lmp_data.r12.data() + size_x12 * 2,
    lmp_data.r12.data() + size_x12 * 3,
    lmp_data.r12.data() + size_x12 * 4,
    lmp_data.r12.data() + size_x12 * 5);
  CUDA_CHECK_KERNEL

  find_descriptor_small_box<<<grid_size, BLOCK_SIZE>>>(
    paramb,
    annmb,
    N,
    N1,
    nlocal,
    nep_data.NN_radial.data(),
    nep_data.NL_radial.data(),
    nep_data.NN_angular.data(),
    nep_data.NL_angular.data(),
    lmp_data.type.data(),
    lmp_data.r12.data(),
    lmp_data.r12.data() + size_x12,
    lmp_data.r12.data() + size_x12 * 2,
    lmp_data.r12.data() + size_x12 * 3,
    lmp_data.r12.data() + size_x12 * 4,
    lmp_data.r12.data() + size_x12 * 5,
    false,//is_polarizability
#ifdef USE_TABLE
    nep_data.gn_radial.data(),
    nep_data.gn_angular.data(),
#endif
    nep_data.potential_per_atom.data(),
    nep_data.Fp.data(),
    nep_data.virial_per_atom.data(),
    nep_data.sum_fxyz.data());
  CUDA_CHECK_KERNEL
  // // bool is_dipole = paramb.model_type == 1;
  find_force_radial_small_box<<<grid_size, BLOCK_SIZE>>>(
    paramb,
    annmb,
    N,
    N1,
    nlocal,
    nep_data.NN_radial.data(),
    nep_data.NL_radial.data(),
    lmp_data.type.data(),
    lmp_data.r12.data(),
    lmp_data.r12.data() + size_x12,
    lmp_data.r12.data() + size_x12 * 2,
    nep_data.Fp.data(),
    false, //is_dipole,
#ifdef USE_TABLE
    nep_data.gnp_radial.data(),
#endif
    nep_data.force_per_atom.data(),
    nep_data.force_per_atom.data() + n_all,
    nep_data.force_per_atom.data() + n_all * 2,
    nep_data.virial_per_atom.data(),
    nep_data.total_virial.data());
  CUDA_CHECK_KERNEL
  find_force_angular_small_box<<<grid_size, BLOCK_SIZE>>>(
    paramb,
    annmb,
    N,
    N1,
    nlocal,
    nep_data.NN_angular.data(),
    nep_data.NL_angular.data(),
    lmp_data.type.data(),
    lmp_data.r12.data() + size_x12 * 3,
    lmp_data.r12.data() + size_x12 * 4,
    lmp_data.r12.data() + size_x12 * 5,
    nep_data.Fp.data(),
    nep_data.sum_fxyz.data(),
    false,//is_dipole,
#ifdef USE_TABLE
    nep_data.gn_angular.data(),
    nep_data.gnp_angular.data(),
#endif
    nep_data.force_per_atom.data(),
    nep_data.force_per_atom.data() + n_all,
    nep_data.force_per_atom.data() + n_all * 2,
    nep_data.virial_per_atom.data(),
    nep_data.total_virial.data());
  CUDA_CHECK_KERNEL
  // checkMemoryUsage(4);
  cudaDeviceSynchronize();

  nep_data.total_virial.copy_to_host(cpu_total_virial);
  nep_data.potential_per_atom.copy_to_host(cpu_potential_per_atom);
  nep_data.force_per_atom.copy_to_host(cpu_force_per_atom);
  // checkMemoryUsage(5);
}


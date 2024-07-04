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

#pragma once
#include "../utilities/common.cuh"
#include "../utilities/gpu_vector.cuh"
#include <tuple>
#include <utility> // for std::move

#define PARAM_SIZE 100

struct LMP_Data  {
  GPU_Vector<float> r12;
  GPU_Vector<int> type;
  GPU_Vector<int> ilist;
  GPU_Vector<int> numneigh;
  GPU_Vector<int> firstneigh;
  GPU_Vector<double> position;
};

struct NEP3_Data {
  GPU_Vector<float> f12x; // 3-body or manybody partial forces
  GPU_Vector<float> f12y; // 3-body or manybody partial forces
  GPU_Vector<float> f12z; // 3-body or manybody partial forces
  GPU_Vector<float> Fp;
  GPU_Vector<float> sum_fxyz;
  GPU_Vector<int> NN_radial;    // radial neighbor list
  GPU_Vector<int> NL_radial;    // radial neighbor list
  GPU_Vector<int> NN_angular;   // angular neighbor list
  GPU_Vector<int> NL_angular;   // angular neighbor list
  GPU_Vector<float> parameters; // parameters to be optimized
  GPU_Vector<int> cell_count;
  GPU_Vector<int> cell_count_sum;
  GPU_Vector<int> cell_contents;

  GPU_Vector<double> potential_per_atom;
  GPU_Vector<double> force_per_atom;
  GPU_Vector<double> virial_per_atom;
  GPU_Vector<double> total_virial;

  std::vector<int> cpu_NN_radial;
  std::vector<int> cpu_NN_angular;
#ifdef USE_TABLE
  GPU_Vector<float> gn_radial;   // tabulated gn_radial functions
  GPU_Vector<float> gnp_radial;  // tabulated gnp_radial functions
  GPU_Vector<float> gn_angular;  // tabulated gn_angular functions
  GPU_Vector<float> gnp_angular; // tabulated gnp_angular functions
#endif
};

class NEP3
{
public:
  struct ParaMB {
    int version = 2; // NEP version, 2 for NEP2 and 3 for NEP3
    int model_type = 0; // 0=potential, 1=dipole, 2=polarizability, 3=temperature-dependent free energy
    float rc_radial = 0.0f;     // radial cutoff
    float rc_angular = 0.0f;    // angular cutoff
    float rcinv_radial = 0.0f;  // inverse of the radial cutoff
    float rcinv_angular = 0.0f; // inverse of the angular cutoff
    int MN_radial = 200;
    int MN_angular = 100;
    int n_max_radial = 0;  // n_radial = 0, 1, 2, ..., n_max_radial
    int n_max_angular = 0; // n_angular = 0, 1, 2, ..., n_max_angular
    int L_max = 0;         // l = 0, 1, 2, ..., L_max
    int dim_angular;
    int num_L;
    int basis_size_radial = 8;  // for nep3
    int basis_size_angular = 8; // for nep3
    int num_types_sq = 0;       // for nep3
    int num_c_radial = 0;       // for nep3
    int num_types = 0;
    float q_scaler[140];
  };

  struct ANN {
    int dim = 0;                 // dimension of the descriptor
    int num_neurons1 = 0;        // number of neurons in the 1st hidden layer
    int num_para = 0;            // number of parameters
    int num_c2 = 0;
    int num_c3 = 0;
    const float* w0[PARAM_SIZE]; // weight from the input layer to the hidden layer
    const float* b0[PARAM_SIZE]; // bias for the hidden layer
    const float* w1[PARAM_SIZE]; // weight from the hidden layer to the output layer
    const float* b1;             // bias for the output layer
    const float* c;
    // for the scalar part of polarizability
    const float* w0_pol[10];
    const float* b0_pol[10];
    const float* w1_pol[10];
    const float* b1_pol;
  };
  NEP3();
  void init_from_file(const char* file_potential, const bool is_rank_0, const int in_device_id);

  ~NEP3(void);
  // int inum;
  // int nlocal;
  // int partition; // max(inum, nlocal)
  // int nghost;
  // int nall;

  ParaMB paramb;
  ANN annmb;
  NEP3_Data nep_data;
  LMP_Data lmp_data;
  std::vector<int> map_atom_type_idx;
  std::vector<int> element_atomic_number_list;
  int atom_nums = 0;
  int atom_nlocal = 0;
  int atom_nums_all = 0;

  void update_potential(float* parameters, ANN& ann);
  // void update_potential_from_cpu(std::vector<float> parameters, ANN& ann);
  void rest_nep_data(int max_atom_nums, int n_local, int n_all, int max_neighbor, bool build_neighbor);

  void checkMemoryUsage(int sgin=0);
#ifdef USE_TABLE
  void construct_table(float* parameters);
#endif
  // void set_partition(int n_all, int n_local, int n_ghost, int n_inum);

  void compute_small_box(
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
    );

  void compute_small_box_optim(
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
    );
  bool has_dftd3 = false;
  bool rank_0 = false;
  int device_id;
};

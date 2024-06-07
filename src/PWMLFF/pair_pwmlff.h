/* -*- c++ -*- ----------------------------------------------------------
     PWmat-MLFF to LAMMPS
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pwmlff, PairPWMLFF);
// clang-format on
#else



#ifndef LMP_PAIR_MLFF_H
#define LMP_PAIR_MLFF_H
#include "nep.h"
#include "pair.h"
#include <iostream>
#include <torch/script.h>
#include <torch/torch.h>


namespace LAMMPS_NS {

    class PairPWMLFF : public Pair {
        public:
            PairPWMLFF(class LAMMPS *);
            ~PairPWMLFF() override;

            int nmax;
            double*** f_n;
            double** e_atom_n;

            std::tuple<std::vector<int>, std::vector<int>, std::vector<double>> generate_neighdata();
            std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<double>> generate_neighdata_nep();
            void compute(int, int) override;
            void settings(int, char **) override;
            void coeff(int, char **) override;
            double init_one(int, int) override;
            void init_style() override;
            int pack_reverse_comm(int, int, double* ) override;
            void unpack_reverse_comm(int, int*, double* ) override;
            void grow_memory();
            std::pair<double, double> calc_max_error(double***, double**);

        protected:
            virtual void allocate();
        
        private:
            int me;
            int num_ff;
            int p_ff_idx;
            unsigned seed;

            NEP3 nep_model;
            std::vector<NEP3> nep_models;

            torch::jit::script::Module module;//dp and nep jit model
            std::vector<torch::jit::script::Module> modules;

            std::vector<double> max_err_list;
            std::vector<double> max_err_ei_list;
            std::string explrError_fname = "explr.error";
            std::FILE *explrError_fp;
            int out_freq = 1;

            torch::Device device = torch::kCPU;
            torch::Dtype dtype = torch::kFloat32;
            std::vector<int> atom_types;
            std::vector<int> model_atom_type_idx;
            int model_ntypes;
            int model_type; // 0 for jitmodel(dp or nep) 1 for nep
            // DP params
            double cutoff;
            // NEP params
            double cutoff_radial;
            double cutoff_angular;
            //common params
            int max_neighbor;
            std::string model_name;

            // std::vector<int> imagetype, imagetype_map, neighbor_list;
            std::vector<int> imagetype_map, neighbor_list;
            std::vector<int> neighbor_type_list; // for nep find neigh and forward
            std::vector<double> dR_neigh;
            // std::vector<int> use_type;

    };

}
#endif
#endif

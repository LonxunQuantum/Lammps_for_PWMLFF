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

#include "pair.h"
#include <iostream>
#include <torch/script.h>
#include <torch/torch.h>


namespace LAMMPS_NS {

    class PairPWMLFF : public Pair {
        public:
            PairPWMLFF(class LAMMPS *);
            ~PairPWMLFF() override;

            std::tuple<std::vector<int>, std::vector<int>, std::vector<double>> generate_neighdata();
            void compute(int, int) override;
            void settings(int, char **) override;
            void coeff(int, char **) override;
            double init_one(int, int) override;
            void init_style() override;
            double calc_max_f_error(std::vector<torch::Tensor>);
            int pack_reverse_comm(int, int, double* ) override;
            void unpack_reverse_comm(int, int*, double* ) override;

        protected:
            virtual void allocate();
        
        private:
            int me;
            int num_ff;
            int p_ff_idx;
            unsigned seed;

            torch::jit::script::Module module;
            std::vector<torch::jit::script::Module> modules;
            std::vector<torch::Tensor> all_forces;
            std::vector<double> max_err_list;
            std::string explrError_fname;
            std::FILE *explrError_fp;
            torch::Device device = torch::kCPU;
            torch::Dtype dtype = torch::kFloat32;
            std::vector<int> type_map;
            double cutoff;
            int max_neighbor;
            std::string model_name;

            // std::vector<int> imagetype, imagetype_map, neighbor_list;
            std::vector<int> imagetype_map, neighbor_list;
            std::vector<double> dR_neigh;
            std::vector<int> use_type;

    };

}
#endif
#endif

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

#include "nep3.cuh"
#include "utilities/common.cuh"
#include "utilities/nep_utilities.cuh"

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 600)
static __device__ __inline__ double atomicAdd(double* address, double val)
{
  unsigned long long* address_as_ull = (unsigned long long*)address;
  unsigned long long old = *address_as_ull, assumed;
  do {
    assumed = old;
    old =
      atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));

  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

static __global__ void find_neighbor_list_small_box(
  NEP3::ParaMB paramb,
  const int n_all,
  const int N, // inum
  const int N1, // 0
  const int N2, // inum
  const int NM,
  const int* __restrict__ g_type,
  const int* __restrict__ g_ilist,
  const int* __restrict__ g_numneigh,
  const int* __restrict__ g_firstneigh,
  const double* __restrict__ g_x,
  const double* __restrict__ g_y,
  const double* __restrict__ g_z,
  int* g_NN_radial,
  int* g_NL_radial,
  int* g_NN_angular,
  int* g_NL_angular,
  float* g_x12_radial,
  float* g_y12_radial,
  float* g_z12_radial,
  float* g_x12_angular,
  float* g_y12_angular,
  float* g_z12_angular)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x + N1;
  if (n1 < N) {
    int count_radial = 0;
    int count_angular = 0;
    int i = g_ilist[n1];
    int jnum = g_numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int n2 = g_firstneigh[n1*NM + jj];
      double x12 = g_x[n1] - g_x[n2];
      double y12 = g_y[n1] - g_y[n2];
      double z12 = g_z[n1] - g_z[n2];
      float distance_square = float(x12 * x12 + y12 * y12 + z12 * z12);

      if (distance_square < paramb.rc_radial * paramb.rc_radial) {
        g_NL_radial[ count_radial * N + n1] = n2;
        g_x12_radial[count_radial * N + n1] = float(x12);
        g_y12_radial[count_radial * N + n1] = float(y12);
        g_z12_radial[count_radial * N + n1] = float(z12);
        count_radial++;
      }
      if (distance_square < paramb.rc_angular * paramb.rc_angular) {
        g_NL_angular[count_angular * N + n1] = n2;
        g_x12_angular[count_angular * N + n1] = float(x12);
        g_y12_angular[count_angular * N + n1] = float(y12);
        g_z12_angular[count_angular * N + n1] = float(z12);
        count_angular++;
      }
    } //for
    g_NN_radial[n1] = count_radial;
    g_NN_angular[n1] = count_angular;
  }// if
}

static __global__ void find_descriptor_small_box(
  NEP3::ParaMB paramb,
  NEP3::ANN annmb,
  const int N,
  const int N1,
  const int nlocal,
  const int* g_NN_radial,
  const int* g_NL_radial,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* __restrict__ g_type,
  const float* __restrict__ g_x12_radial,
  const float* __restrict__ g_y12_radial,
  const float* __restrict__ g_z12_radial,
  const float* __restrict__ g_x12_angular,
  const float* __restrict__ g_y12_angular,
  const float* __restrict__ g_z12_angular,
  const bool is_polarizability,
#ifdef USE_TABLE
  const float* __restrict__ g_gn_radial,
  const float* __restrict__ g_gn_angular,
#endif
  double* g_pe,
  float* g_Fp,
  double* g_virial,
  float* g_sum_fxyz)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x + N1;
  if (n1 < N) {
    int t1 = g_type[n1];
    float q[MAX_DIM] = {0.0f};
    // get radial descriptors
    for (int i1 = 0; i1 < g_NN_radial[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL_radial[index];
      float r12[3] = {g_x12_radial[index], g_y12_radial[index], g_z12_radial[index]};
      float d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
#ifdef USE_TABLE
      int index_left, index_right;
      float weight_left, weight_right;
      find_index_and_weight(
        d12 * paramb.rcinv_radial, index_left, index_right, weight_left, weight_right);
      int t12 = t1 * paramb.num_types + g_type[n2];
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        q[n] +=
          g_gn_radial[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_left +
          g_gn_radial[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_right;
      }
#else
      float fc12;
      find_fc(paramb.rc_radial, paramb.rcinv_radial, d12, fc12);
      int t2 = g_type[n2];
      float fn12[MAX_NUM_N];
      if (paramb.version == 2) {
        find_fn(paramb.n_max_radial, paramb.rcinv_radial, d12, fc12, fn12);
        for (int n = 0; n <= paramb.n_max_radial; ++n) {
          float c = (paramb.num_types == 1)
                      ? 1.0f
                      : annmb.c[(n * paramb.num_types + t1) * paramb.num_types + t2];
          q[n] += fn12[n] * c;
        }
      } else {
        find_fn(paramb.basis_size_radial, paramb.rcinv_radial, d12, fc12, fn12);
        for (int n = 0; n <= paramb.n_max_radial; ++n) {
          float gn12 = 0.0f;
          for (int k = 0; k <= paramb.basis_size_radial; ++k) {
            int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
            c_index += t1 * paramb.num_types + t2;
            gn12 += fn12[k] * annmb.c[c_index];
          }
          q[n] += gn12;
        }
      }
#endif
    }

    // get angular descriptors
    for (int n = 0; n <= paramb.n_max_angular; ++n) {
      float s[NUM_OF_ABC] = {0.0f};
      for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
        int index = i1 * N + n1;
        int n2 = g_NL_angular[index];
        float r12[3] = {g_x12_angular[index], g_y12_angular[index], g_z12_angular[index]};
        float d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
#ifdef USE_TABLE
        int index_left, index_right;
        float weight_left, weight_right;
        find_index_and_weight(
          d12 * paramb.rcinv_angular, index_left, index_right, weight_left, weight_right);
        int t12 = t1 * paramb.num_types + g_type[n2];
        float gn12 =
          g_gn_angular[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n] *
            weight_left +
          g_gn_angular[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n] *
            weight_right;
        accumulate_s(d12, r12[0], r12[1], r12[2], gn12, s);
#else
        float fc12;
        find_fc(paramb.rc_angular, paramb.rcinv_angular, d12, fc12);
        int t2 = g_type[n2];
        if (paramb.version == 2) {
          float fn;
          find_fn(n, paramb.rcinv_angular, d12, fc12, fn);
          fn *=
            (paramb.num_types == 1)
              ? 1.0f
              : annmb.c
                  [((paramb.n_max_radial + 1 + n) * paramb.num_types + t1) * paramb.num_types + t2];
          accumulate_s(d12, r12[0], r12[1], r12[2], fn, s);
        } else {
          float fn12[MAX_NUM_N];
          find_fn(paramb.basis_size_angular, paramb.rcinv_angular, d12, fc12, fn12);
          float gn12 = 0.0f;
          for (int k = 0; k <= paramb.basis_size_angular; ++k) {
            int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
            c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
            gn12 += fn12[k] * annmb.c[c_index];
          }
          accumulate_s(d12, r12[0], r12[1], r12[2], gn12, s);
        }
#endif
      }
      if (paramb.num_L == paramb.L_max) {
        find_q(paramb.n_max_angular + 1, n, s, q + (paramb.n_max_radial + 1));
      } else if (paramb.num_L == paramb.L_max + 1) {
        find_q_with_4body(paramb.n_max_angular + 1, n, s, q + (paramb.n_max_radial + 1));
      } else {
        find_q_with_5body(paramb.n_max_angular + 1, n, s, q + (paramb.n_max_radial + 1));
      }
      for (int abc = 0; abc < NUM_OF_ABC; ++abc) {
        g_sum_fxyz[(n * NUM_OF_ABC + abc) * nlocal + n1] = s[abc];
      }
    }

    // nomalize descriptor
    for (int d = 0; d < annmb.dim; ++d) {
      q[d] = q[d] * paramb.q_scaler[d];
    }

    // get energy and energy gradient
    float F = 0.0f, Fp[MAX_DIM] = {0.0f};

    if (is_polarizability) {
      apply_ann_one_layer(
        annmb.dim,
        annmb.num_neurons1,
        annmb.w0_pol[t1],
        annmb.b0_pol[t1],
        annmb.w1_pol[t1],
        annmb.b1_pol,
        q,
        F,
        Fp,
        t1);
      // Add the potential values to the diagonal of the virial
      g_virial[n1] = F;
      g_virial[n1 + N * 1] = F;
      g_virial[n1 + N * 2] = F;

      F = 0.0f;
      for (int d = 0; d < annmb.dim; ++d) {
        Fp[d] = 0.0f;
      }
    }

    apply_ann_one_layer(
      annmb.dim, annmb.num_neurons1, annmb.w0[t1], annmb.b0[t1], annmb.w1[t1], annmb.b1, q, F, Fp, t1);
    g_pe[n1] += F;
    for (int d = 0; d < annmb.dim; ++d) {
      g_Fp[d * nlocal + n1] = Fp[d] * paramb.q_scaler[d];
    }
    // if (n1 == N-1){
    //   printf("g_pe[%d] = %f g_Fp_d0[%d] = %f\n", n1, g_pe[n1], n1, g_Fp[0 * nlocal + n1]);
    // }
  }
}

static __global__ void find_force_radial_small_box(
  NEP3::ParaMB paramb,
  NEP3::ANN annmb,
  const int N,
  const int N1,
  const int nlocal,
  const int* g_NN,
  const int* g_NL,
  const int* __restrict__ g_type,
  const float* __restrict__ g_x12,
  const float* __restrict__ g_y12,
  const float* __restrict__ g_z12,
  const float* __restrict__ g_Fp,
  const bool is_dipole,
#ifdef USE_TABLE
  const float* __restrict__ g_gnp_radial,
#endif
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_total_virial)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x + N1;
  if (n1 < N) {
    // if (n1 == 0){
    //   printf("find_force_radial_small_box n1 = 0\n");
    // }
    int t1 = g_type[n1];
    for (int i1 = 0; i1 < g_NN[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL[index];
      int t2 = g_type[n2];
      float r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      float d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      float d12inv = 1.0f / d12;
      float f12[3] = {0.0f};
#ifdef USE_TABLE
      int index_left, index_right;
      float weight_left, weight_right;
      find_index_and_weight(
        d12 * paramb.rcinv_radial, index_left, index_right, weight_left, weight_right);
      int t12 = t1 * paramb.num_types + t2;
      for (int n = 0; n <= paramb.n_max_radial; ++n) {
        float gnp12 =
          g_gnp_radial[(index_left * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_left +
          g_gnp_radial[(index_right * paramb.num_types_sq + t12) * (paramb.n_max_radial + 1) + n] *
            weight_right;
        float tmp12 = g_Fp[n1 + n * nlocal] * gnp12 * d12inv;
        for (int d = 0; d < 3; ++d) {
          f12[d] += tmp12 * r12[d];
        }
      }
#else
      float fc12, fcp12;
      find_fc_and_fcp(paramb.rc_radial, paramb.rcinv_radial, d12, fc12, fcp12);
      float fn12[MAX_NUM_N];
      float fnp12[MAX_NUM_N];
      if (paramb.version == 2) {
        find_fn_and_fnp(paramb.n_max_radial, paramb.rcinv_radial, d12, fc12, fcp12, fn12, fnp12);
        for (int n = 0; n <= paramb.n_max_radial; ++n) {
          float tmp12 = g_Fp[n1 + n * nlocal] * fnp12[n] * d12inv;
          tmp12 *= (paramb.num_types == 1)
                     ? 1.0f
                     : annmb.c[(n * paramb.num_types + t1) * paramb.num_types + t2];
          for (int d = 0; d < 3; ++d) {
            f12[d] += tmp12 * r12[d];
          }
        }
      } else {
        find_fn_and_fnp(
          paramb.basis_size_radial, paramb.rcinv_radial, d12, fc12, fcp12, fn12, fnp12);
        for (int n = 0; n <= paramb.n_max_radial; ++n) {
          float gnp12 = 0.0f;
          for (int k = 0; k <= paramb.basis_size_radial; ++k) {
            int c_index = (n * (paramb.basis_size_radial + 1) + k) * paramb.num_types_sq;
            c_index += t1 * paramb.num_types + t2;
            gnp12 += fnp12[k] * annmb.c[c_index];
          }
          float tmp12 = g_Fp[n1 + n * nlocal] * gnp12 * d12inv;
          for (int d = 0; d < 3; ++d) {
            f12[d] += tmp12 * r12[d];
          }
        }
      }
#endif
      double s_sxx = 0.0;
      double s_sxy = 0.0;
      double s_sxz = 0.0;
      double s_syx = 0.0;
      double s_syy = 0.0;
      double s_syz = 0.0;
      double s_szx = 0.0;
      double s_szy = 0.0;
      double s_szz = 0.0;
      if (is_dipole) {
        double r12_square = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
        s_sxx -= r12_square * f12[0];
        s_syy -= r12_square * f12[1];
        s_szz -= r12_square * f12[2];
      } else {
        s_sxx -= r12[0] * f12[0];
        s_syy -= r12[1] * f12[1];
        s_szz -= r12[2] * f12[2];
      }
      s_sxy -= r12[0] * f12[1];
      s_sxz -= r12[0] * f12[2];
      s_syz -= r12[1] * f12[2];
      s_syx -= r12[1] * f12[0];
      s_szx -= r12[2] * f12[0];
      s_szy -= r12[2] * f12[1];

      atomicAdd(&g_fx[n1], double(f12[0]));
      atomicAdd(&g_fy[n1], double(f12[1]));
      atomicAdd(&g_fz[n1], double(f12[2]));
      atomicAdd(&g_fx[n2], double(-f12[0]));
      atomicAdd(&g_fy[n2], double(-f12[1]));
      atomicAdd(&g_fz[n2], double(-f12[2]));
      // save virial
      // xx xy xz    0 3 4
      // yx yy yz    6 1 5
      // zx zy zz    7 8 2
      atomicAdd(&g_virial[n2 + 0 * N], s_sxx);
      atomicAdd(&g_virial[n2 + 1 * N], s_syy);
      atomicAdd(&g_virial[n2 + 2 * N], s_szz);
      atomicAdd(&g_virial[n2 + 3 * N], s_sxy);
      atomicAdd(&g_virial[n2 + 4 * N], s_sxz);
      atomicAdd(&g_virial[n2 + 5 * N], s_syz);
      atomicAdd(&g_virial[n2 + 6 * N], s_syx);
      atomicAdd(&g_virial[n2 + 7 * N], s_szx);
      atomicAdd(&g_virial[n2 + 8 * N], s_szy);

      atomicAdd(&g_total_virial[0], s_sxx);// xx
      atomicAdd(&g_total_virial[1], s_syy);// yy
      atomicAdd(&g_total_virial[2], s_szz);// zz
      atomicAdd(&g_total_virial[3], s_sxy);// xy
      atomicAdd(&g_total_virial[4], s_sxz);// xz
      atomicAdd(&g_total_virial[5], s_syz);// yz
    }
  }
}

static __global__ void find_force_angular_small_box(
  NEP3::ParaMB paramb,
  NEP3::ANN annmb,
  const int N,
  const int N1,
  const int nlocal,
  const int* g_NN_angular,
  const int* g_NL_angular,
  const int* __restrict__ g_type,
  const float* __restrict__ g_x12,
  const float* __restrict__ g_y12,
  const float* __restrict__ g_z12,
  const float* __restrict__ g_Fp,
  const float* __restrict__ g_sum_fxyz,
  const bool is_dipole,
#ifdef USE_TABLE
  const float* __restrict__ g_gn_angular,
  const float* __restrict__ g_gnp_angular,
#endif
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_total_virial)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x + N1;
  if (n1 < N) {
    // if (n1 == 0){
    //   printf("find_force_angular_small_box n1 = 0\n");
    // }
    float Fp[MAX_DIM_ANGULAR] = {0.0f};
    float sum_fxyz[NUM_OF_ABC * MAX_NUM_N];
    for (int d = 0; d < paramb.dim_angular; ++d) {
      Fp[d] = g_Fp[(paramb.n_max_radial + 1 + d) * nlocal + n1];
    }
    for (int d = 0; d < (paramb.n_max_angular + 1) * NUM_OF_ABC; ++d) {
      sum_fxyz[d] = g_sum_fxyz[d * nlocal + n1];
    }

    int t1 = g_type[n1];

    for (int i1 = 0; i1 < g_NN_angular[n1]; ++i1) {
      int index = i1 * N + n1;
      int n2 = g_NL_angular[n1 + N * i1];
      float r12[3] = {g_x12[index], g_y12[index], g_z12[index]};
      float d12 = sqrt(r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2]);
      float f12[3] = {0.0f};

#ifdef USE_TABLE
      int index_left, index_right;
      float weight_left, weight_right;
      find_index_and_weight(
        d12 * paramb.rcinv_angular, index_left, index_right, weight_left, weight_right);
      int t12 = t1 * paramb.num_types + g_type[n2];
      for (int n = 0; n <= paramb.n_max_angular; ++n) {
        int index_left_all =
          (index_left * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n;
        int index_right_all =
          (index_right * paramb.num_types_sq + t12) * (paramb.n_max_angular + 1) + n;
        float gn12 =
          g_gn_angular[index_left_all] * weight_left + g_gn_angular[index_right_all] * weight_right;
        float gnp12 = g_gnp_angular[index_left_all] * weight_left +
                      g_gnp_angular[index_right_all] * weight_right;
        if (paramb.num_L == paramb.L_max) {
          accumulate_f12(n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp, sum_fxyz, f12);
        } else if (paramb.num_L == paramb.L_max + 1) {
          accumulate_f12_with_4body(
            n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp, sum_fxyz, f12);
        } else {
          accumulate_f12_with_5body(
            n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp, sum_fxyz, f12);
        }
      }
#else
      float fc12, fcp12;
      find_fc_and_fcp(paramb.rc_angular, paramb.rcinv_angular, d12, fc12, fcp12);
      int t2 = g_type[n2];
      if (paramb.version == 2) {
        for (int n = 0; n <= paramb.n_max_angular; ++n) {
          float fn;
          float fnp;
          find_fn_and_fnp(n, paramb.rcinv_angular, d12, fc12, fcp12, fn, fnp);
          const float c =
            (paramb.num_types == 1)
              ? 1.0f
              : annmb.c
                  [((paramb.n_max_radial + 1 + n) * paramb.num_types + t1) * paramb.num_types + t2];
          fn *= c;
          fnp *= c;
          accumulate_f12(n, paramb.n_max_angular + 1, d12, r12, fn, fnp, Fp, sum_fxyz, f12);
        }
      } else {
        float fn12[MAX_NUM_N];
        float fnp12[MAX_NUM_N];
        find_fn_and_fnp(
          paramb.basis_size_angular, paramb.rcinv_angular, d12, fc12, fcp12, fn12, fnp12);
        for (int n = 0; n <= paramb.n_max_angular; ++n) {
          float gn12 = 0.0f;
          float gnp12 = 0.0f;
          for (int k = 0; k <= paramb.basis_size_angular; ++k) {
            int c_index = (n * (paramb.basis_size_angular + 1) + k) * paramb.num_types_sq;
            c_index += t1 * paramb.num_types + t2 + paramb.num_c_radial;
            gn12 += fn12[k] * annmb.c[c_index];
            gnp12 += fnp12[k] * annmb.c[c_index];
          }
          if (paramb.num_L == paramb.L_max) {
            accumulate_f12(n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp, sum_fxyz, f12);
          } else if (paramb.num_L == paramb.L_max + 1) {
            accumulate_f12_with_4body(
              n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp, sum_fxyz, f12);
          } else {
            accumulate_f12_with_5body(
              n, paramb.n_max_angular + 1, d12, r12, gn12, gnp12, Fp, sum_fxyz, f12);
          }
        }
      }
#endif
      double s_sxx = 0.0;
      double s_sxy = 0.0;
      double s_sxz = 0.0;
      double s_syx = 0.0;
      double s_syy = 0.0;
      double s_syz = 0.0;
      double s_szx = 0.0;
      double s_szy = 0.0;
      double s_szz = 0.0;
      if (is_dipole) {
        double r12_square = r12[0] * r12[0] + r12[1] * r12[1] + r12[2] * r12[2];
        s_sxx -= r12_square * f12[0];
        s_syy -= r12_square * f12[1];
        s_szz -= r12_square * f12[2];
      } else {
        s_sxx -= r12[0] * f12[0];
        s_syy -= r12[1] * f12[1];
        s_szz -= r12[2] * f12[2];
      }
      s_sxy -= r12[0] * f12[1];
      s_sxz -= r12[0] * f12[2];
      s_syz -= r12[1] * f12[2];
      s_syx -= r12[1] * f12[0];
      s_szx -= r12[2] * f12[0];
      s_szy -= r12[2] * f12[1];

      atomicAdd(&g_fx[n1], double(f12[0]));
      atomicAdd(&g_fy[n1], double(f12[1]));
      atomicAdd(&g_fz[n1], double(f12[2]));
      atomicAdd(&g_fx[n2], double(-f12[0]));
      atomicAdd(&g_fy[n2], double(-f12[1]));
      atomicAdd(&g_fz[n2], double(-f12[2]));
      // save virial
      // xx xy xz    0 3 4
      // yx yy yz    6 1 5
      // zx zy zz    7 8 2
      atomicAdd(&g_virial[n2 + 0 * N], s_sxx);
      atomicAdd(&g_virial[n2 + 1 * N], s_syy);
      atomicAdd(&g_virial[n2 + 2 * N], s_szz);
      atomicAdd(&g_virial[n2 + 3 * N], s_sxy);
      atomicAdd(&g_virial[n2 + 4 * N], s_sxz);
      atomicAdd(&g_virial[n2 + 5 * N], s_syz);
      atomicAdd(&g_virial[n2 + 6 * N], s_syx);
      atomicAdd(&g_virial[n2 + 7 * N], s_szx);
      atomicAdd(&g_virial[n2 + 8 * N], s_szy);

      atomicAdd(&g_total_virial[0], s_sxx);// xx
      atomicAdd(&g_total_virial[1], s_syy);// yy
      atomicAdd(&g_total_virial[2], s_szz);// zz
      atomicAdd(&g_total_virial[3], s_sxy);// xy
      atomicAdd(&g_total_virial[4], s_sxz);// xz
      atomicAdd(&g_total_virial[5], s_syz);// yz
    }
  }
}


/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef GRAN_SUB_MOD_CLASS
// clang-format off
GranSubModStyle(none,
         GranSubModNormalNone,
         NORMAL);

GranSubModStyle(hooke,
         GranSubModNormalHooke,
         NORMAL);

GranSubModStyle(hertz,
         GranSubModNormalHertz,
         NORMAL);

GranSubModStyle(hertz/material,
         GranSubModNormalHertzMaterial,
         NORMAL);

GranSubModStyle(dmt,
         GranSubModNormalDMT,
         NORMAL);

GranSubModStyle(jkr,
         GranSubModNormalJKR,
         NORMAL);
// clang-format on
#else

#ifndef GRAN_SUB_MOD_NORMAL_H
#define GRAN_SUB_MOD_NORMAL_H

#include "gran_sub_mod.h"

namespace LAMMPS_NS {
namespace Granular_NS {

class GranSubModNormal : public GranSubMod {
 public:
  GranSubModNormal(class GranularModel *, class LAMMPS *);
  ~GranSubModNormal() {};
  virtual bool touch();
  virtual double pulloff_distance(double, double);
  virtual double calculate_contact_radius();
  virtual double calculate_forces() = 0;
  virtual void set_fncrit();
  double damp; // argument historically needed by damping
  double Emod, poiss;
  double Fncrit;
  int material_properties, cohesive_flag;
};

/* ---------------------------------------------------------------------- */

class GranSubModNormalNone : public GranSubModNormal {
 public:
  GranSubModNormalNone(class GranularModel *, class LAMMPS *);
  double calculate_forces();
};

/* ---------------------------------------------------------------------- */

class GranSubModNormalHooke : public GranSubModNormal {
 public:
  GranSubModNormalHooke(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  double calculate_forces();
 protected:
  double k;
};

/* ---------------------------------------------------------------------- */

class GranSubModNormalHertz : public GranSubModNormal {
 public:
  GranSubModNormalHertz(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  double calculate_forces();
 protected:
  double k;
};

/* ---------------------------------------------------------------------- */

class GranSubModNormalHertzMaterial : public GranSubModNormalHertz {
 public:
  GranSubModNormalHertzMaterial(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void mix_coeffs(double*, double*) override;
};

/* ---------------------------------------------------------------------- */

class GranSubModNormalDMT : public GranSubModNormal {
 public:
  GranSubModNormalDMT(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void mix_coeffs(double*, double*) override;
  double calculate_forces();
  void set_fncrit() override;
 protected:
  double k, cohesion;
  double F_pulloff, Fne;
};

/* ---------------------------------------------------------------------- */

class GranSubModNormalJKR : public GranSubModNormal {
 public:
  GranSubModNormalJKR(class GranularModel *, class LAMMPS *);
  void coeffs_to_local() override;
  void mix_coeffs(double*, double*) override;
  bool touch() override;
  double pulloff_distance(double, double);
  double calculate_contact_radius() override;
  double calculate_forces();
  void set_fncrit() override;
 protected:
  double k, cohesion;
  double Emix, F_pulloff, Fne;
};

}    // namespace Granular_NS
}    // namespace LAMMPS_NS

#endif /*GRAN_SUB_MOD_NORMAL_H */
#endif /*GRAN_SUB_MOD_CLASS_H */

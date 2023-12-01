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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(dipole,AtomVecDipole);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_DIPOLE_H
#define LMP_ATOM_VEC_DIPOLE_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecDipole : virtual public AtomVec {
 public:
  AtomVecDipole(class LAMMPS *);

  void grow_pointers() override;
  void data_atom_post(int) override;

 protected:
  double **mu;
};

}    // namespace LAMMPS_NS

#endif
#endif

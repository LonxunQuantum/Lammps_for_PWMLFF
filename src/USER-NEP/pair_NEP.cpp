/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Junjie Wang (Nanjing University)
------------------------------------------------------------------------- */

#include "pair_NEP.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "nep.h"
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>

#define LAMMPS_VERSION_NUMBER 20230208 // use the new neighbor list starting from this version

using namespace LAMMPS_NS;

PairNEP::PairNEP(LAMMPS* lmp) : Pair(lmp)
{
#if LAMMPS_VERSION_NUMBER >= 20201130
  centroidstressflag = CENTROID_AVAIL;
#else
  centroidstressflag = 2;
#endif

  restartinfo = 0;
  manybody_flag = 1;

  single_enable = 0;

  inited = false;
  allocated = 0;
}

PairNEP::~PairNEP()
{
  if (copymode)
    return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

void PairNEP::allocate()
{
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      setflag[i][j] = 1;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  allocated = 1;
}

void PairNEP::coeff(int narg, char** arg)
{
  if (!allocated)
    allocate();
}

void PairNEP::settings(int narg, char** arg)
{
  if (narg != 1) {
    error->all(FLERR, "Illegal pair_style command; nep requires 1 parameter");
  }
  model_filename = arg[0];
}

void PairNEP::init_style()
{
#if LAMMPS_VERSION_NUMBER >= 20220324
  neighbor->add_request(this, NeighConst::REQ_FULL);
#else
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
#endif

  bool is_rank_0 = (comm->me == 0);
  nep_model.init_from_file(model_filename, is_rank_0);
  inited = true;
  cutoff = nep_model.paramb.rc_radial;
  cutoffsq = cutoff * cutoff;
  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
      cutsq[i][j] = cutoffsq;
}

double PairNEP::init_one(int i, int j) { return cutoff; }

void PairNEP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag, vflag);
  }
  double total_potential = 0.0;
  double total_virial[6] = {0.0};
  double* per_atom_potential = nullptr;
  double** per_atom_virial = nullptr;
  if (eflag_atom) {
    per_atom_potential = eatom;
  }
  if (cvflag_atom) {
    per_atom_virial = cvatom;
  }

  nep_model.compute_for_lammps(
   atom->nlocal, list->inum, list->ilist, list->numneigh, list->firstneigh, atom->type, atom->x, total_potential,
    total_virial, per_atom_potential, atom->f, per_atom_virial, 0);

  // int nghost = atom->nghost;
  // for (int i = 0; i < list->inum + nghost; i++)
  // {
  //   std::cout <<"=== force " << i << " " << atom->f[i][0] << " " << atom->f[i][1] << " " << atom->f[i][2] << std::endl;
  // }

  if (eflag) {
    eng_vdwl += total_potential;
    // std::cout << "=== total " << total_potential << std::endl;
  }
  // std::cout << "=== total_potential " << total_potential << std::endl;
  if (vflag) {
    // std::cout<< "=== total_virial 6 " << total_virial[0] << " " << total_virial[1] << " "  << total_virial[2] << " "  << total_virial[3] << " "  << total_virial[4] << " "  << total_virial[5] << std::endl; 
    for (int component = 0; component < 6; ++component) {
      virial[component] += total_virial[component];
    }
  //  std::cout<< "=== virial 6 " << virial[0] << " " << virial[1] << " "  << virial[2] << " "  << virial[3] << " "  << virial[4] << " "  << virial[5] << std::endl; 
  }
}

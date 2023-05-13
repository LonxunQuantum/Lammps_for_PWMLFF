/*
    Interface to PWMATMLFF
    2023.1
*/

#include "pair_pwmatmlff.h"

//#include "diy.h"   

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
//#include "respa.h"
#include "update.h"
#include "domain.h"

#include <cmath>
#include <cstring>
// tmp
#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <string>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std; 

// declration of f2c! 

extern "C" 
{
    void f2c_calc_energy_force( int* /*imodel*/,
                                int* /*nlocal*/,
                                int* /*nlocal + nghost*/, 
                                int* /*ntypes*/,
                                int* /*type*/,
                                int* /*itype*/, 
                                int* /*num_neigh*/,
                                int* /*list_neigh*/,
                                double* /*dR_neigh*/,
                                double* /*e_tot*/, 
                                double* /*f*/, 
                                double* /*virial*/);
    // ff data loading 
    void dp_ff_load(char*, int*, int*);

    // ff data destroy  
    // pointer of the name string 
    void dp_ff_deallocate(int* ); 
}
    
PairPWMATMLFF::PairPWMATMLFF(LAMMPS *lmp) : Pair(lmp)
{
    writedata = 1;
}


PairPWMATMLFF::~PairPWMATMLFF()
{
    if (copymode) return;

    /*
        dump stat file
    */
    string name = "explr.stat";
    ofstream file;

    file.open(name, ios::trunc);
    file << num_success << " ";
    file << num_cand << " ";
    file << num_fail << " "; 
    file.close(); 
    
    /*
        dump error list 
    */
    string err_name = "explr.error";
    file.open(err_name);

    for(int i=0; i< max_err_list.size(); i++)
        file << max_err_list[i] <<endl; 
    file.close(); 

    //cout << "num of candidate:" << atom_config_list.size() <<endl; 
    /*
        pick out 10 config for later selection
    */

    //int num_cand = atom_config_list.size(); 
    int num_select = 10;  

    if (atom_config_list.size() > 0 )
    {
        vector<int> idx(num_cand,0); 
        
        for (int i=0; i<num_cand; i++)
            idx[i] = i; 
        
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

        shuffle(idx.begin(),idx.end(), std::default_random_engine(seed)); 
        
        for(int i=0; i< num_select; i++ )
        {
            /*
                write candidates' configs
                RANDOMLY CHOOSE 10 CONFIGS
                do padding if num_cand < num_select    
            */
            if (i<num_cand)
                atom_config_list[idx[i]].write(i);   
            else
                atom_config_list[idx[num_cand-1]].write(i);   
        }
    }

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(cut);
        memory->destroy(cutsq);
        memory->destroy(scale);
        memory->destroy(e_atom);
        memory->destroy(itype_atom);
        memory->destroy(num_neigh);
        memory->destroy(list_neigh);
        memory->destroy(dR_neigh);
    }

    // ff paras memory 
    for (int ff_idx = 1; ff_idx <= num_ff; ff_idx ++)
    {
        cout << "!!!Releasing force field memory " << ff_idx << "!!!"<<endl;
        dp_ff_deallocate(&ff_idx);
        cout << "Force field memory released" << endl;
    }
    
}   

void PairPWMATMLFF::allocate()
{   
    allocated = 1;  
    int n = atom->ntypes + 1;

    memory->create(setflag, n, n, "pair_pwmatmlff:setflag");

    for (int i = 1; i < n; i++)
        for (int j = i; j < n; j++) 
            setflag[i][j] = 0;
    
    memory->create(cut,n,n,"pair_pwmatmlff:cut");
    memory->create(cutsq,n,n,"pair_pwmatmlff:cutsq");
    memory->create(scale,n,n,"pair_pwmatmlff:scale");
    
    memory->create(e_atom, atom->natoms, "pair_pwmatmlff:e_atom");
    memory->create(itype_atom, atom->natoms, "pair_pwmatmlff:itype_atom");

    memory->create(num_neigh, atom->natoms, atom->ntypes, "pair_pwmatmlff:num_neigh");  
    memory->create(list_neigh, atom->natoms, atom->ntypes, 100, "pair_pwmatmlff:list_neigh");
    memory->create(dR_neigh, atom->natoms, atom->ntypes, 100, 3, "pair_pwmatmlff:dR_neigh");
    
    atom_config_list.reserve(500);

}

double PairPWMATMLFF::calc_max_f_err(int nlocal, double** f, double** f_n)
{
    /*
        f : force of reference model 
        f_n: force of all other models 
    */
    vector<double> f_ave;           
    vector<double> f_error[num_ff];

    f_ave.resize(3*nlocal);
        
    // start from model #2 
    for(int i =0; i<num_ff; i++)
        f_error[i].resize(3*nlocal);

    // copy from model #1
    for(int i=0; i<nlocal; i++)
    {
        f_ave[3*i] = f[i][0];
        f_ave[3*i+1] = f[i][1];
        f_ave[3*i+2] = f[i][2]; 
    }
    //cout << "copy done" <<endl; 

    // sum over all models
    for (int ff_idx = 1; ff_idx < num_ff; ff_idx++)
    {
        for (int j = 0; j < 3*nlocal; j++)
            f_ave[j] += f_n[ff_idx][j];
    }

    // calc ensemble average
    double q = 1.0/num_ff; 
    for (int i=0; i< 3*nlocal; i++)
        f_ave[i] *= q; 
    
    //cout << "average done" <<endl; 

    // calc error 
        
    for (int j=0; j < nlocal; j++)
    {
        f_error[0][j*3] = f[j][0] - f_ave[j*3];
        f_error[0][j*3+1] = f[j][1] - f_ave[j*3+1];
        f_error[0][j*3+2] = f[j][2] - f_ave[j*3+2];
    }

    for (int idx = 1; idx< num_ff; idx++)
        for(int j=0; j<3*nlocal; j++)
        {
            f_error[idx][j] = f_n[idx][j] - f_ave[j]; 
        }
        
    double max_err = -1.0; 

    //cout << "error done" <<endl; 
        
    // find max error 
    for (int idx = 0; idx< num_ff; idx++)
        for(int j=0; j<3*nlocal; j+=3)
        {
            double err = f_error[idx][j]  *  f_error[idx][j] 
                       + f_error[idx][j+1] * f_error[idx][j+1]    
                       + f_error[idx][j+2] * f_error[idx][j+2];

            err = sqrt(err); 

            if (err > max_err) 
                max_err = err; 
        } 
        
    //cout << "max error : " << max_err <<endl; 
    return max_err; 
}

void PairPWMATMLFF::compute(int eflag, int vflag)
{

  double* h = domain->h;
  double** x = atom->x;
  double** f = atom->f;
  int* type = atom->type;
  int* tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int ntypes = atom->ntypes;
  double* virial = force->pair->virial;

  int debug = 0;

  // for benchmark in exploration 
  double e_tot_tmp; 
  double virial_tmp[6]; 
  
  double* f_n[num_ff];
  double* e_atom_n[num_ff];   

  for(int i=1; i<num_ff; i++)
  {
      f_n[i] = new double[3*nlocal];
      e_atom_n[i] = new double[nlocal];
  }    

  ev_init(eflag, vflag);

  // get lattice
  lattice[0] = h[0];   // xx
  lattice[4] = h[1];   // yy
  lattice[8] = h[2];   // zz
  lattice[7] = h[3];   // yz
  lattice[6] = h[4];   // xz
  lattice[3] = h[5];   // xy

  for (int i = 0; i < nlocal; i++)
  {
    // atom->type starts from 1 in lammps.data,
    // and type_map index starts from 0;
    itype_atom[i] = type_map[type[i]-1];
  }

  generate_neighdata();
  
  // call fortran subroutine energy force 
  /*
      do calc_energy_force for 4 models 
      1 is used for MD
      2, 3, 4 for the test
  */  
  
  for (int ff_idx = 1; ff_idx <= num_ff; ff_idx++)
  {   
      // usual MD 
      if (ff_idx == 1)
      {
          f2c_calc_energy_force(  &imodel, 
                                  &nlocal,
                                  &nall, 
                                  &ntypes,
                                  &itype_atom[0], 
                                  &num_neigh[0][0],
                                  &list_neigh[0][0][0],
                                  &dR_neigh[0][0][0][0],
                                  &e_atom[0], 
                                  &e_tot,
                                  &f[0][0], 
                                  &virial[0],
                                  &ff_idx);
      // benchmark of the rest models 
      } else {
          f2c_calc_energy_force(  &imodel, 
                                  &nlocal,
                                  &nall, 
                                  &ntypes,
                                  &itype_atom[0], 
                                  &num_neigh[0][0],
                                  &list_neigh[0][0][0],
                                  &dR_neigh[0][0][0][0],
                                  &e_atom_n[ff_idx-1], 
                                  &e_tot_tmp,
                                  &f[ff_idx-1], 
                                  &virial_tmp[0],
                                  &ff_idx);
      }
  }


  if (debug == 1) {
    // this reverse comm is only for test
    //comm->reverse_comm();
    printf("@@@ ----- f2c end --------- @@@\n");
    printf("@@@ force begin @@@\n\n");
    for (int i = 0; i < nlocal; i++)
    {
      //printf("@@@ %7d %7d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n", i, tag[i], \
      //            f[i][0], f[i][1], f[i][2], x[i][0], x[i][1], x[i][2]);
      printf("@@@ %7d %7d %12.6f %12.6f %12.6f\n", i, tag[i], \
                  f[i][0], f[i][1], f[i][2]);

    }
    //printf("@@@ energy begin @@@\n\n");
    //for (int i = 0; i < nlocal; i++)
    //{
    //  printf("@@@ %7d %7d %12.6f\n", i, tag[i], e_atom[i]);
    //}
  }

  //printf("@@@ virial %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", virial[0],virial[1],virial[2],virial[3],virial[4],virial[5]);

  /*
      exploration mode.
      select candidates
      write as fractional coordinate 
      Note: only write data at the very end (in the destructor)! 
  */  
  if (num_ff != 1)
  {   
      // calculate model deviation with Force 
      double max_err = calc_max_f_err(nlocal, f, f_n); 
      
      max_err_list.push_back(max_err);

      //cout << "max error:" << max_err <<endl; 
      
      if (max_err > candidate_err &&  max_err < fail_err)
      //if (true)
      {   
          num_cand ++; 

          atom_config tmp(nlocal);
          tmp.get_all(lattice, x, itype_atom); 
          //tmp.print_all(); 
          atom_config_list.push_back(tmp);
          //save image for SCF; 
      } 
      else if (max_err > fail_err)
      {
          num_fail++; 
      }
      else if (max_err < candidate_err)
      {
          num_success++; 
      }
  }   

  // deallocate benchmark arrays 
  for(int i=1; i<num_ff; i++)
  {
      delete [] f_n[i];
      delete [] e_atom_n[i];
  } 

  // potential energy
  if (eflag_global)
    eng_vdwl += e_tot;
        
}

void PairPWMATMLFF::settings(int argc, char** argv)
{
}

void PairPWMATMLFF::coeff(int narg, char** arg)
{
    if (!allocated) allocate();

    int ntypes = atom->ntypes;

    //if (narg != (2+2+ntypes*2+4+2)) {    
    //    error->one(FLERR, "pair_coeff error:\
    //      * * imodel elem1 elem2 ... elem_ntypes cut1 cut2 ... cut_ntypes");
    //}
    
    imodel = utils::inumeric(FLERR, arg[2], false, lmp);
    num_ff = utils::inumeric(FLERR, arg[3], false, lmp);

    for (int i = 0; i < ntypes; i++)
      type_map[i] = utils::inumeric(FLERR, arg[i+2+2+num_ff], false, lmp);
   
    for (int i = 1; i <= ntypes; i++) {
      setflag[i][i] = 1;
      cut[i][i] = utils::numeric(FLERR, arg[i+2+2+num_ff+ntypes], false, lmp);
    }

    for (int ff_idx=1; ff_dix < 5; ff_idx++) {

      int name_len=0; 
      for (int i=0; i< 500; i++) {
        if (arg[2+2+ff_idx-1][i] == '\0')
          break; 
        name_len++; 
      }

      if (comm->me == 0)
        printf("\nStart loading force field: \n\n%s %d\n\n", \
          arg[narg-1], name_len);

      dp_ff_load(arg[2+2+ff_idx-1], &ff_idx, &name_len);

      if (comm->me == 0)
        printf("Force field loaded successfully\n\n"); 
    }

    if (narg > (2+2+num_ff+ntypes*2)) {
      candidate_err = utils::numeric(FLERR, arg[narg-2], false, lmp);
      fail_err = utils::numeric(FLERR, arg[narg-1], false, lmp);
    }

    num_fail = 0;
    num_cand = 0;
    num_success = 0;
}

double  PairPWMATMLFF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }
  return cut[i][j];
}

void PairPWMATMLFF::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR, "Pair style PWMATMLFF requires newon pair on");

  neighbor->add_request(this, NeighConst::REQ_FULL);
}


void PairPWMATMLFF::generate_neighdata()
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double** x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int etnum;
  for (i = 0; i < nlocal; i++)
    for (j = 0; j < ntypes; j++)
      num_neigh[i][j] = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    //printf("i %7d %7d %12.6f\n", i, itype, cutsq[itype][itype]);
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = x[j][0] - xtmp;
      dely = x[j][1] - ytmp;
      delz = x[j][2] - ztmp;
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][itype]) {
        etnum = num_neigh[i][jtype-1];
        dR_neigh[i][jtype-1][etnum][0] = delx;
        dR_neigh[i][jtype-1][etnum][1] = dely;
        dR_neigh[i][jtype-1][etnum][2] = delz;
        // fortran start from 1
        list_neigh[i][jtype-1][etnum] = j+1;
        num_neigh[i][jtype-1] += 1;
        //printf("  j %7d %7d %7d \n", j, jtype, num_neigh[i][jtype-1]);
      }
    }
  }

  /*
  for (i = 0; i < nlocal; i++) {
    printf("i   %7d\n", i);
    for (int t = 0; t < ntypes; t++) {
      printf("  t   %7d %7d\n", t, num_neigh[i][t]);
      for (j = 0; j < num_neigh[i][t]; j++)
        printf("      %9.6f %9.6f %9.6f \n", \
          dR_neigh[i][t][j][0], dR_neigh[i][t][j][1], dR_neigh[i][t][j][2]);
    }
  }
  */

}

//void Pair::init_list(int /*which*/, NeighList *ptr)
//{ 
//  list = ptr;
//}




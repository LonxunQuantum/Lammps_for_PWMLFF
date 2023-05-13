/* -*- c++ -*- ----------------------------------------------------------
    liuliping, PWmat-MLFF to LAMMPS
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(pwmatmlff,PairPWMATMLFF);
// clang-format on
#else

#ifndef LMP_PAIR_PWMATMLFF_H
#define LMP_PAIR_PWMATMLFF_H

#include "pair.h"
#include <vector>
#include <unordered_map>

using namespace std;

typedef unsigned long long ull; 

namespace LAMMPS_NS {

class atom_config
{
    /*
        a single atom.config 
    */
public:
    int num_atom; 
    vector<double> lattice; 
    vector<double> position;
    vector<int> type;
    
    atom_config(int num)
    {
        num_atom = num; 
        lattice.resize(9);
        position.resize(3*num_atom);
        type.resize(num_atom);
    }
    
    void get_all(double* latt, double** x , int* itype_atom)
    {
        for (int i=0; i< 9; i++)
            lattice[i] = latt[i];
        
        for (int i=0; i< num_atom; i++)
        {
            position[i*3] = x[i][0];
            position[i*3+1] = x[i][1];
            position[i*3+2] = x[i][2]; 

            type[i] = itype_atom[i]; 
        }
        
    }
    void write(int idx)
    {
        /*
            write int atom.config style 
        */
        string name = "config."+to_string(idx);
        
        ofstream file;
        
        file.precision(12); 
        
        file.open(name);
        
        file << num_atom<<endl; 
        file << "LATTICE" <<endl; 

        for (int i=0; i<9; i++)
        {
            if (i%3 ==0 && i!=0)
                file << endl; 
            file << setw(22) << lattice[i];
        }    
        
        file << endl << "POSITION" <<endl; 
        
        for (int i =0; i < num_atom; i++)
        {
            file << setw(3) << type[i]; 
            file << setw(22) << position[3*i] << " ";
            file << setw(22) << position[3*i+1] << " ";
            file << setw(22) << position[3*i+2]; 
            file << "       0      0       0" <<endl ;  
        }

        file.close();
    }   

    void print_all()
    {
        cout.precision(12); 

        cout << num_atom <<endl; 

        for (int i=0; i< 9; i++) 
        {
            if (i%3 ==0 )
                cout << endl; 
            cout << lattice[i] << " ";
        }
        cout << endl; 
        
        for (int i =0; i < num_atom; i++)
        {
            cout << type[i] << " "; 
            cout << position[3*i] << " ";
            cout << position[3*i+1] << " ";
            cout << position[3*i+2] << endl;
        }
    }
};


class PairPWMATMLFF : public Pair {

public:
    /*
        1) arrays that will be passed into fortran are all consequetive chunks.  
           get_fort_idx_xd functions take care of the indices, so that elements are accessed
           in the fortran style

        2) use STL containers smartly. Less stress on memory management

        3) 
    */
    
    PairPWMATMLFF(class LAMMPS *);
    ~PairPWMATMLFF() override;
    
    int num_atoms;
    int ievery = 1000;
    int iago = ievery;
    int imodel = 1; // should be read from model file
    int iflag_reneighbor = 1;
    double** scale;
    double** cut;
    
    int* itype_atom;
    int* itype_tmp;
    double* e_atom;
    double lattice[9];
    double reclat[9];
    double tmp_v3[3];
    double e_tot;
    int** num_neigh;
    int*** list_neigh;
    double**** dR_neigh;

    // read type_map from pair_coeff of lmp.in
    // max 10 elements, this array should be [29, 8, 0,...] for CuO
    
    int type_map[10];
    int num_ff;                      // number of force field 

    int num_fail;
    int num_cand;
    int num_success; 

    double fail_err = 0.25;
    double candidate_err = 0.05; 
    
    vector<double> max_err_list; 
    vector<atom_config> atom_config_list; 
    
    // major member functions 
    double calc_max_f_err(int, double**, double**);
    void compute(int, int) override;
    void settings(int, char **) override;
    void coeff(int, char **) override;
    void init_style() override;
    double init_one(int, int) override;
    void read_ff(string); 
    virtual void allocate();
    void generate_neighdata();
    
};

}    // namespace LAMMPS_NS

#endif
#endif

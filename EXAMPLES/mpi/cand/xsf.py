
import os
from  units import mass_list, element_index, charge_list
from read_write_conf import \
       Read_One_Conf, \
       Write_One_PWmat_Data, \
       Write_One_LAMMPS_Data, \
       Cart2Lamda, Lamda2cart 
import numpy as np
import copy
from tri2orth import Tri2Orth

def Composition(elem):

  natoms = len(elem)
  atype = np.zeros(natoms, dtype=int)
  charge = np.zeros(natoms, dtype=int)
  eindex = np.zeros(natoms, dtype=int)
  elem_list = list(dict.fromkeys(elem))
  num_atype = len(elem_list)  
  mass = np.zeros(num_atype, dtype=int)

  for i in range(natoms):
    eindex[i] = element_index[elem[i]]
    atype[i] = elem_list.index(elem[i]) + 1
    charge[i] = charge_list[elem[i]]

  for i in range(num_atype):
    mass[i] = mass_list[elem_list[i]]

  return num_atype, atype, mass, charge, eindex, elem_list

ddir = "./"

for cid in range(1000):

  filename = ddir + "/cand_" + str(cid) + ".config"

  infile = open(filename, 'r')
  lines = infile.readlines()

  print("cid %7d" % cid)

  lattice, Etot, natoms, elem, coords, forces, endln = Read_One_Conf(lines)

  infile.close()

  num_atype, atype, mass, charge, eindex, elem_list = Composition(elem)
  #coords = Cart2Lamda(lattice, coords)
    
  # ----------------------------------------------------
  #  write out pwmat movement
  # ----------------------------------------------------
  filename = './0' + str(cid)
  sofile = open(filename, 'w')
  
  Write_One_PWmat_Data(sofile,cid,natoms,Etot,lattice,eindex,atype,charge,coords,forces)
  sofile.close()
  
  # ----------------------------------------------------
  #  write out lammps data
  # ----------------------------------------------------
  filename = './' + str(cid) + '.data'
  
  #lattice[0], lattice[1], lattice[2] = Tri2Orth(lattice[0], lattice[1], lattice[2])
  coords = Lamda2cart(lattice, coords)
  Write_One_LAMMPS_Data(filename, natoms, lattice, num_atype, mass, atype, coords, charge)




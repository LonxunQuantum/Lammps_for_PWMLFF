
import numpy as np
from units import index_element

# -------------------------------------------------
#  change real coordinates to fraction coordinates
# -------------------------------------------------
def Cart2Lamda(lattice, coords):

  matrix = np.array(lattice).reshape((3,3))
  inv_matrix = np.linalg.inv(matrix)

  frac_coord = []
  for r in range(len(coords)):
    frac_coord.append(np.dot(coords[r],inv_matrix))

  return frac_coord

# -------------------------------------------------
#  change fraction coordinates to real coordinates
# -------------------------------------------------
def Lamda2cart(lattice, xs):

  xhi = lattice[0][0]
  yhi = lattice[1][1]
  zhi = lattice[2][2]
  xy  = lattice[1][0]
  xz  = lattice[2][0]
  yz  = lattice[2][1]

  Cart = []
  for i in range(len(xs)):
    tx = xs[i][0]*xhi + xs[i][1]*xy  + xs[i][2]*xz
    ty =                xs[i][1]*yhi + xs[i][2]*yz
    tz =                               xs[i][2]*zhi
    Cart.append([tx,ty,tz])

  return Cart

# -------------------------------------------------
#  read one PWmat final config formation configuration
# -------------------------------------------------
def Read_One_Conf(lines):

  lattice = []
  Etot = 0.0
  natoms = 0
  elem = []
  coords = []
  forces = []

  natoms = int(lines[0].split()[0])

  endln = 0
  for i, line in enumerate(lines):

    if "Lattice" in line:
      for j in range(i+1, i+1+3):
        lattice.append([float(c) for c in lines[j].split()])

    if "Position" in line:
      for j in range(i+1, i+1+natoms):
        tokens = lines[j].split()
        Z = int(tokens[0])
        elem.append(index_element[Z])
        coords.append([float(j) for j in tokens[1:4]])

      forces = np.zeros([natoms,3], dtype=float)

      return lattice, Etot, natoms, elem, coords, forces, endln


# -------------------------------------------------
#  write one lammps data formation configuration
# -------------------------------------------------
def Write_One_LAMMPS_Data(filename, natoms, lattice, num_atype, \
  mass, atype, coords, charge=[], action='w'):

  xhi = lattice[0][0]
  yhi = lattice[1][1]
  zhi = lattice[2][2]
  xy  = lattice[1][0]
  xz  = lattice[2][0]
  yz  = lattice[2][1]
  xlo = ylo = zlo = 0.0

  ofile = open(filename, action)

  ofile.write("only for test\n\n")
  ofile.write("%9d atoms\n" % (natoms))
  ofile.write("%5d atom types\n\n" % (num_atype))
  ofile.write("%16.12f %16.12f xlo xhi\n" % (xlo, xhi))
  ofile.write("%16.12f %16.12f ylo yhi\n" % (ylo, yhi))
  ofile.write("%16.12f %16.12f zlo zhi\n" % (zlo, zhi))
  ofile.write("%16.12f %16.12f %16.12f xy xz yz\n\n" % (xy, xz, yz))

  ofile.write("Masses\n\n")

  for i in range(num_atype):
    ofile.write("%7d  %16.9f\n" % (i+1, mass[i]))

  conftype = 1
  if len(charge) == natoms:
    conftype = 2
  if conftype == 1:
    ofile.write("\nAtoms # atomic\n\n")
    for i in range(natoms):
      ofile.write("%9d %5d %16.12f %16.12f %16.12f\n" \
        % (i+1, atype[i], coords[i][0], coords[i][1], coords[i][2]))
  elif conftype == 2:
    ofile.write("\nAtoms # charge\n\n")
    for i in range(natoms):
      ofile.write("%9d %5d %5d %16.12f %16.12f %16.12f\n" \
        % (i+1, atype[i], charge[i], coords[i][0], coords[i][1], coords[i][2]))

  ofile.close()


# -------------------------------------------------
#  write one pwmat data formation configuration
# -------------------------------------------------
def Write_One_PWmat_Data(sofile,conf_id,natoms,Etot,lattice,\
  elem,atype,charge,coords,f):

  sofile.write("%9d atoms,Iteration (fs) = %9d, Etot,Ep,Ek (eV) = %16.9f %16.9f\n" %\
    (natoms, int(conf_id), Etot, Etot))

  sofile.write("Lattice vector (Angstrom)\n")
  for l in range(3):
    sofile.write("%12.9f %12.9f %12.9f\n" % \
      (lattice[l][0], lattice[l][1], lattice[l][2]))

  sofile.write("Position (nomalized)\n")
  for i in range(natoms):  
    sofile.write("%7d %7d %5d %12.9f %12.9f %12.9f 1 1 1\n" % \
      (elem[i], atype[i], charge[i], coords[i][0], coords[i][1], coords[i][2]))

  sofile.write("Force (-force, eV/Angstrom)\n")
  for i in range(natoms):
    sofile.write("%7d %12.9f %12.9f %12.9f\n" % \
      (elem[i], f[i][0], f[i][1], f[i][2]))

  sofile.write("Atomic-Energy, Etot(ev), E_nonloc(eV),Q_atom:dE(eV)= %12.6f\n" % Etot)
  per_atom_e = Etot/natoms
  for i in range(natoms):
    sofile.write("%7d %12.9f %12.9f %12.9f\n" % \
      (elem[i], per_atom_e, 0.0, 0.0))

  sofile.write("-------------------------------------------------\n")


# -------------------------------------------------
#  write one MOVEMENT formation configuration
# -------------------------------------------------
def Write_One_MOVEMENT(sofile,conf_id,natoms,Etot,lattice,\
  elem,atype,charge,coords,f):

  sofile.write("%9d atoms,Iteration (fs) = %9d, Etot,Ep,Ek (eV) = %16.9f %16.9f\n" %\
    (natoms, int(conf_id), Etot, Etot))

  sofile.write("Lattice vector (Angstrom)\n")
  for l in range(3):
    sofile.write("%12.9f %12.9f %12.9f\n" % \
      (lattice[l][0], lattice[l][1], lattice[l][2]))

  sofile.write("Position (nomalized)\n")
  for i in range(natoms):  
    sofile.write("%7d %12.9f %12.9f %12.9f 1 1 1\n" % \
      (elem[i], coords[i][0], coords[i][1], coords[i][2]))

  sofile.write("Force (-force, eV/Angstrom)\n")
  for i in range(natoms):
    sofile.write("%7d %12.9f %12.9f %12.9f\n" % \
      (elem[i], f[i][0], f[i][1], f[i][2]))

  sofile.write("Atomic-Energy, Etot(ev), E_nonloc(eV),Q_atom:dE(eV)= %12.6f\n" % Etot)
  per_atom_e = Etot/natoms
  for i in range(natoms):
    sofile.write("%7d %12.9f %12.9f %12.9f\n" % \
      (elem[i], per_atom_e, 0.0, 0.0))

  sofile.write("-------------------------------------------------\n")




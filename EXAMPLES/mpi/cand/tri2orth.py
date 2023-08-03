
import numpy as np

def normalize(v):
  norm = np.linalg.norm(v)
  if norm == 0:
    return [0, 0, 0]
  return v/norm

## general lattice vector
# A --> | Ax Ay Az |
# B --> | Bx By Bz |
# C --> | Cx Cy Cz |
## ensure the lattice vector a along x direction
#                            b in the x-y plane
#                            a, b, c form the right hand axis
# a --> | ax  0  0 |      | xx  0  0 |
# b --> | bx by  0 | <==> | xy yy  0 |
# c --> | cx cy cz |      | xz yz zz |
def Tri2Orth(A,B,C):

  nA = normalize(A)
  Al = np.linalg.norm(A)
  Bl = np.linalg.norm(B)
  Cl = np.linalg.norm(C)

  ax = np.linalg.norm(A)
  bx = np.dot(B,nA)
  by = np.sqrt(Bl*Bl-bx*bx)
  cx = np.dot(C,nA)
  cy = (np.dot(B,C)-bx*cx)/by
  cz = np.sqrt(Cl*Cl - cx*cx - cy*cy)

  xx = ax
  yy = by
  zz = cz
  xy = bx
  xz = cx
  yz = cy

  rx = [xx, 0.0, 0.0]
  ry = [xy,  yy, 0.0]
  rz = [xz,  yz, zz]

  a, b, c = rx, ry, rz 
  
  xlo = 0.0
  xhi = a[0]
  ylo = 0.0
  yhi = b[1]
  zlo = 0.0
  zhi = c[2]
  xy  = b[0]
  xz  = c[0]
  yz  = c[1]

  #return xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
  return a, b, c 



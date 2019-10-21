import numpy as np
import glob

from tqdm import tqdm

from ase.io import read,write
from ase import Atom,Atoms

def unique_check(material,rc):
  cell = material.get_cell()
  kk = []
  for m in range(len(material)):
	for n in range(m+1,len(material)):
	  dmn = material.get_distance(m,n,mic=True)
	  if dmn < rc:
		kk.append(n)
  kk = np.array(kk); new_atom = []
  for m in range(len(material)):
	det = (kk==m)
	if not det.any():
	  new_atom.append(material[m])

  new_mat = Atoms(new_atom,cell=cell,pbc=[1,1,1])
  return new_mat


for cucumber in tqdm(glob.glob('*.vasp')):
  org_mat = read(cucumber)
  eles = org_mat.get_chemical_symbols()
  mat_O = org_mat.copy(); mat_V = org_mat.copy()
  (a,) = np.where(np.array(eles)=='O'); del mat_V[a]
  (b,) = np.where(np.array(eles)=='V'); del mat_O[b]

  mat_V2 = unique_check(mat_V,1.6)
  mat_O2 = unique_check(mat_O,1.3)

  for mm in mat_O2:
	mat_V2.append(mm)

  write(cucumber,mat_V2)

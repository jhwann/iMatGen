import glob
import numpy as np
from tqdm import tqdm
from ase.io import read,write
from ase import Atoms,Atom

def partial_check(atoms):
  k = []; N = len(atoms)
  for idx_i in range(N):
	k.append([])
	for idx_j in range(idx_i+1,N):
	  dij = atoms.get_distance(idx_i,idx_j,mic=True)
	  if dij < 1.2:
		k[-1].append(idx_j)

  return k

def atoms_expansion(atoms,k):
  new_atoms = atoms.copy()
  ll = [];
  for idx,k_ele in enumerate(k):
	if len(k_ele) > 0:
	  if atoms[k_ele[0]].symbol == 'O':
		ll.append(k_ele[0])
	  elif atoms[idx] == 'O':
		ll.append(idx)

	ll = list(set(ll))
  del new_atoms[ll]
  return new_atoms


def neighbors(atoms,idx,k_ele):
  dic = {'V':0,'O':0}
  lists = range(len(atoms))
  lists.remove(idx); lists.remove(k_ele)
  for ii in lists:
	dij = atoms.get_distance(idx,ii,mic=True)
	if dij <= 2.0:
	  dic[atoms[ii].symbol] += 1

  return dic

for name in tqdm(glob.glob('*.vasp')):
	atoms = read(name)
	partial_list = np.array(partial_check(atoms))
	new_atoms = atoms_expansion(atoms,partial_list)
	while True:
		partial_list = np.array(partial_check(atoms))
		new_atoms = atoms_expansion(atoms,partial_list)
		if atoms.get_chemical_formula()==new_atoms.get_chemical_formula():
			break
	write(name,new_atoms)

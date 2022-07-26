import numpy as np

import ase
from ase.io import read, write
from ase import Atoms
#from ase.build import molecule

#from ase.calculators.emt import EMT
from gpaw import GPAW, PW

from ase.optimize import BFGS

#filename = 'from_ase_database.xyz'
filename = 'addH_reorder_opt_test2_indices.xyz'

# read xyz file
mol = read(filename, format = 'xyz')
#mol = molecule('CH3CH2OH')
mol.center(vacuum=6.0)

print(mol.get_positions())
print(mol.cell)

# add calculator
#calc = EMT()
calc = GPAW(xc='PBE', mode=PW(400), txt='out.txt')
mol.calc = calc

# compute potential energy
e = mol.get_potential_energy()
print('Pontential energy:', e)

# optimize molecule
opt = BFGS(mol, trajectory = 'opt.traj', logfile='opt.log')
opt.run(fmax = 0.01)

e_opt = mol.get_potential_energy()
print('Pontential energy after optimization:', e_opt)

# write structure
write('opt_' + filename, mol)

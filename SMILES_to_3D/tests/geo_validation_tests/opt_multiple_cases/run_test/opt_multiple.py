import os
import numpy as np
import time

import ase
from ase.io import read, write
from ase import Atoms
from gpaw import GPAW, PW
from ase.optimize import BFGS

cwd = os.getcwd()
output = open('summary.txt', 'w')
output.write('ID    energy before opt    energy after opt \n')

for i in range(10):
    print('Molecule:', i)
    root_folder = os.path.join(cwd, str(i))
    in_filename = os.path.join(root_folder, str(i) + '.xyz')

    # read xyz file
    mol = read(in_filename, format = 'xyz')
    # add unit cell, distance between two repeated molecule is vacuum
    mol.center(vacuum=4.0)
    mol.set_pbc((False, False, False))

    # add calculator
    out_filename = os.path.join(root_folder, str(i) + '_out.txt')
    calc = GPAW(xc='PBE', mode=PW(400), txt=out_filename)
    mol.calc = calc

    # compute potential energy before opt
    e = mol.get_potential_energy()

    # optimize molecule
    traj_filename = os.path.join(root_folder, str(i) + '_opt.traj')
    log_filename = os.path.join(root_folder, str(i) + '_opt.log')
    opt = BFGS(mol, trajectory = traj_filename, logfile=log_filename)
    opt.run(fmax = 0.05, steps=1000)
 
    # compute potential energy after opt
    e_opt = mol.get_potential_energy()

    # write structure
    opt_filename = os.path.join(root_folder, str(i) + '_opt.xyz')
    write(opt_filename, mol)

    # write to summary
    outstr = '{}    {:6.6f}    {:6.6f} \n'.format(i, e, e_opt)
    output.write(outstr)

    # wait
    time.sleep(20)

output.close()

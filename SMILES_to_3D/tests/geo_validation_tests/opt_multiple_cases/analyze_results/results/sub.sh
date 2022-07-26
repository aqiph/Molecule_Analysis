#!/bin/sh
#$ -S /bin/bash
#$ -q lenovo
#$ -pe lenovo 32
#$ -cwd
#$ -e error
#$ -o outpt
#$ -N test

module load mpi/openmpi3-x86_64
MPIRUN=/usr/lib64/openmpi3/bin/mpirun

$MPIRUN -np $NSLOTS python3 example.py



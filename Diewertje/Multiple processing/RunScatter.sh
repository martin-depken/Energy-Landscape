#!/bin/bash
#PBS -l nodes=4:ppn=2

module load mpi/openmpi-1.8.8-gnu
mpirun -np 8 python /home/dddekker/BEP/scatter_MPI.py






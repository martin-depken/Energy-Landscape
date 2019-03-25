#!/bin/bash
#PBS -l nodes=2:ppn=3

module load mpi/openmpi-1.8.8-gnu
/opt/ud/openmpi-1.8.8/bin/mpirun -np 6 python /home/dddekker/BEP/Hello_World_MPI.py






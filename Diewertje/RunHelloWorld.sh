#!/bin/bash
#PBS -l nodes=3:ppn=3

numbr=$((PBS_NUM_PPN*PBS_NUM_NODES))

module load mpi/openmpi-1.8.8-gnu
/home/dddekker/miniconda3/envs/dev/bin/mpirun -np $numbr -hostfile $PBS_NODEFILE python /home/dddekker/BEP/Hello_World_MPI.py








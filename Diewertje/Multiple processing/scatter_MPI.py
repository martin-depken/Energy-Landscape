# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 13:37:40 2019

@author: Diewertje
"""
from mpi4py import MPI

name = MPI.Get_processor_name()

comm = MPI.COMM_WORLD
size=comm.Get_size()
rank=comm.Get_rank()
#print(name[0])

if rank == 0: 
   data = [(x+1) **2 for x in range (size)]
   print('scattering data',data)
else:
   data = None
data = comm.scatter(data,root=0)
print('rank',rank,'has data: ', data, 'on node', name)


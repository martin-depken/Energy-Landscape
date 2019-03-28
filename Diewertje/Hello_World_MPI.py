# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 09:39:33 2019

@author: Diewertje
"""

from mpi4py import MPI

comm = MPI.COMM_WORLD

size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

print ("Hello world, I am rank ", rank," of ", size, " running processes on node ",name,"\n")


# =============================================================================
# Opslaan in de file werkt niet omdat elke core een file aanmaakt, dus ik krijg hem alleen te zien van de core die hem het laatst terug stuurt
# path='/home/dddekker/BEP/Hello_output.txt'
# f=open(path, 'w',1)
# print ("Hello world, I am rank ", rank," of ", size, " running processes on node ",name,"\n", file=f)
# f.close()
# =============================================================================



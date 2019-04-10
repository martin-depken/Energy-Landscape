# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 09:30:12 2019

@author: Diewertje
"""
from scipy import *
from pylab import *
import numpy as np
import Calculate_ABA_Finkelsteinlab_Diewertje as CalcABA
import matplotlib.pylab as plt
import plotting_Boyle_Diewertje as plt_B
import CRISPR_free_energy_landscape_Diewertje as free

#f=open("parameters_Boyle.txt","r")
#parameters=f.read().splitlines()
#parameters=list(map(float,parameters))
filename='../data/22_3_2019/fit_22_3_2019_sim_13.txt'
model_id = 'general_energies_no_kPR'
# filename = simset[np.argsort(chi_squared)[3]]
Nparams = 43
parameters = plt_B.load_simm_anneal(filename, Nparams)

epsilon = parameters[:-2] # for commanted file above it is -3
guide_length=20
F=[]
for i in range(1,guide_length+1): # PAM is 0, that needs to match, so start at 1
    mismatch_positions=np.arange(i,21)
    energies=CalcABA.get_energies(epsilon,mismatch_positions, guide_length)
    F.append(-np.log(np.sum(np.exp(-np.cumsum(energies)))))
    print(mismatch_positions)

#mismatch_positions=[]
#energies=CalcABA.get_energies(epsilon,mismatch_positions, guide_length)
#F_ontarget=-np.log(np.sum(np.exp(-np.cumsum(energies))))
#F=np.array(F)-F_ontarget

plt.plot(range(0,20),F,marker='o') # plot from 0, since you start the mm_block at 1, so that is the free energy of the PAM
plt.xlabel('Start of mismatch block')
plt.ylabel('Free energy')
# plot does not go till 20, since that would be perfect match

_,freeE=free.plot_free_energy_landscape(parameters,model_id,show_plot=False)
plt.plot(range(0,21),freeE,marker='o')

plt.xticks(range(1,21))
plt.grid()
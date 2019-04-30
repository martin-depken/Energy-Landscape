# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:05:16 2019

@author: Diewertje
"""
from scipy import *
from pylab import *
import numpy as np

import sys
sys.path.append('../code_general/') # added this after we run it on the cluster
from read_model_ID import unpack_parameters

# this is for the assumption that we only have 2 states, sol and PAM
def calc_Pbnd(parameters, concentrations, reference,mismatch_positions, model_id = 'general_energies', guide_length = 20, T=10*60):
    # T: evaluate solution at time T
    # Pbnd=vector with the Pbnd at the different concnetrations
    # concentrations= a vector of lenth 8
    # Kon= vector of length 8
    # Koff= 1 value
    rel_concentrations = concentrations/reference
    epsilon, forward_rates = unpack_parameters(parameters, model_id, guide_length)
    #epsilonPAM=epsilon[0]
    epsilonPAM=epsilon[0]+np.log(10) # since parameters are at 
    epsilonCI=np.array([get_energies(epsilon,mismatch_positions, guide_length=20)]) #-1*np.cumsum(epsilon[1:21])+epsilon[0]
    epsilonLandscape=np.cumsum(epsilonCI)
    Ksp0=forward_rates[0]/10 # rate solution to PAM, /10 because paramters are at 10 nM
    Kon=rel_concentrations*Ksp0
    print(Kon)
    Kps=Ksp0*np.exp(epsilonPAM) 
    Ppam=np.exp(-epsilonPAM)/sum(np.exp(-epsilonLandscape))
    Koff=Kps*Ppam
    print(Koff)
    print(type(Koff))
    print(Kon+Koff)
    Pbnd=(Kon/(Kon+Koff))*(1-np.exp(-(Kon+Koff)*T)) 
    return Pbnd
    
def get_energies(epsilon,mismatch_positions, guide_length=20):
    '''
    For general (position dependent) model make a list with the energies at every bound state
    At positions with a mismatch incorporated: add mismatch penalty (epsI[mm_pos])

    So this function returns the minima in the energy lanscape (the actual energy at every state)

    :param epsilon: [epsPAM, epsC[state_1],....,epsC[state_N],epsI[state_1],...epsI[state_N] ]
    provide as np.array()
    :param mismatch_positions: each mismach position has a range [1,2, ... , 20]
    :return: vector with the minima in the landscape
    '''
    if type(mismatch_positions)==type([]):
        mismatch_positions = np.array(mismatch_positions)
    new_epsilon = epsilon.copy()
    epsI = np.array(new_epsilon[(guide_length+1):])
    energies = -1*np.array(new_epsilon[0:(guide_length+1)]) # convention: epsC>0 means downward slope
    energies[0] = new_epsilon[0]                 # convention: epsPAM>0 means upward slope
    if len(mismatch_positions)>0:
        energies[mismatch_positions.astype(int)] += epsI[(mismatch_positions.astype(int)-1)]
    energies=np.ndarray.tolist(energies)
    return energies
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:22:35 2019

@author: Diewertje
"""
from scipy import *
from pylab import *
import numpy as np 
import pandas as pd
import Bio
from Bio.pairwise2 import format_alignment
import copy

#import sys 
#sys.path.append('../code_general_Finkelsteinlab/')

#import Process_SeqLibrary_Finkelsteinlab as preprocess

#import plotting_Finkelsteinlab as pltData

#%matplotlib inline 
import matplotlib.pylab as plt 

import Prepare_data as Pre
import Chisq_Finkelstein as Chi
import SimulatedAnnealing_Finkelstein_parallel as SA
import Calculate_ABA_Finkelsteinlab_Diewertje as ABA
import functools

def main():
    path = '../Data_ABA_Finkelsteinlab/' 
    filename= 'cas9-target-e-replicate-1-delta-abas_Canonical_OT-r_0-2.csv'
    xdata,ydata,yerr=Pre.Prepare_Cdata(path,filename)
    
    fit_result_file= 'fit_results.txt'
    monitor_file = 'monitor.txt'
    init_monitor_file='init_monitor.txt'
    
    model_ID =  'init_limit_general_energies_v2' #'general_energies' FOR GENERAL ENERGIES MAG LAATSTE LOWERBOUND NIET NEGATIEF ZIJN, DUS ZET OP 0
    concentrations = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300]) # in nanoMolair
    reference=1 # in nanomolair
    
    upper_bnd =  [10.0] + [10.0]*40 +  [300.0] *3 #*2
    lower_bnd = [0.0] + [-10.0]*40 + [-7.0] *3 #*2 
    initial_guess = [5.0] + [1.0]*40 + [0.5] *3 #np.loadtxt('parameters.txt')
    
    KineticModel = functools.partial(Chi.calc_Chi_square,model_id=model_ID, guide_length=20,
                                     concentrations=concentrations,reference=reference)
    
    ONtarget=functools.partial(ABA.calc_ABA,concentrations=concentrations,reference=reference, 
                               mismatch_positions=[],model_id=model_ID, guide_length = 20, T=10*60)
    
    print('now starts with the sim_anneal funciton')
    fit_result = SA.sim_anneal_fit(xdata=xdata[:2],
                                    ydata=ydata[:2],
                                    yerr = yerr[:2],
                                    Xstart= np.array(initial_guess),
                                    lwrbnd= np.array(lower_bnd),
                                    upbnd= np.array(upper_bnd),
                                    model='I_am_using_multi_processing_in_stead',
                                    objective_function=KineticModel,
                                    on_target_function=ONtarget,
                                    Tstart=100.,             # infered from run on my computer/other runs on cluster
                                    use_relative_steps=False,
                                    delta=1.0,
                                    tol=1E-5,
                                    Tfinal=50,
                                    adjust_factor=1.1,
                                    cooling_rate=0.6,
                                    N_int=100,
                                    AR_low=40,
                                    AR_high=60,
                                    use_multiprocessing=True,
                                    nprocs=1,
                                    output_file_results = fit_result_file,
                                    output_file_monitor = monitor_file,
                                    output_file_init_monitor=init_monitor_file
                                  )
    
if __name__ == "__main__": 
    main()
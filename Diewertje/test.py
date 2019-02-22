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
import sys
sys.path.append('../code_ABA_Finkelsteinlab/')
import Calculate_ABA_Finkelsteinlab as ABA

import functools

def main():
    path = '../Data_ABA_Finkelsteinlab/' 
    filename= 'cas9-target-e-replicate-1-delta-abas_Canonical_OT-r_0-2.csv'
    xdata,ydata,yerr=Pre.Prepare_Cdata(path,filename)
    
    fit_result_file= 'fit_results.txt'
    monitor_file = 'monitor.txt'
    init_monitor_file='init_monitor.txt'
    
    gRNA_length = 20
    model_ID = 'general_energies'
    parameters =  np.array([5.0] + [1.0]*40 + [0.5] *2)
    concentrations = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300]) # in nanoMolair
    reference=1 # in nanomolair
    
    upper_bnd =  [10.0] + [10.0]*40 +  [3.0] *2
    lower_bnd = [0.0] + [-10.0]*40 + [-7.0] *2
    initial_guess =  [5.0] + [0.0]*40 + [0.0] *2
    
    KineticModel = functools.partial(Chi.calc_Chi_square,model_id=model_ID, guide_length=20,
                                     concentrations=concentrations,reference=reference)
    print('now starts with the sim_anneal funciton')
    fit_result = SA.sim_anneal_fit(xdata=xdata[:2],
                                    ydata=ydata[:2],
                                    yerr = yerr[:2],
                                    Xstart= np.array(initial_guess),
                                    lwrbnd= np.array(lower_bnd),
                                    upbnd= np.array(upper_bnd),
                                    model='I_am_using_multi_processing_in_stead',
                                    objective_function=KineticModel,
                                    on_target_function=ABA.calc_ABA,
                                    Tstart=100.,             # infered from run on my computer/other runs on cluster
                                    use_relative_steps=False,
                                    delta=1.0,
                                    tol=1E-5,
                                    Tfinal=50,
                                    adjust_factor=1.1,
                                    cooling_rate=0.5,
                                    N_int=10,
                                    AR_low=0,
                                    AR_high=100,
                                    use_multiprocessing=True,
                                    nprocs=1,
                                    output_file_results = fit_result_file,
                                    output_file_monitor = monitor_file,
                                    output_file_init_monitor=init_monitor_file
                                  )
    
if __name__ == "__main__": 
    main()
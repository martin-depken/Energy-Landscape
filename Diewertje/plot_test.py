# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 15:21:37 2019

@author: Diewertje
"""
from scipy import *
from pylab import *

import importlib as imp
import matplotlib.pylab as plt 
import numpy as np 
import pandas as pd 
import seaborn as sns
sns.set_style('ticks');
current_colors = sns.color_palette()
sns.set_palette('Accent');

# Adjust all text and axes colors to black. default is some studip gray that appears to faint when printed 
from matplotlib import rc , rcParams
rc('axes',edgecolor='black')
rc('xtick',color='black',labelsize=15)
rc('ytick',color='black',labelsize=15)
rc('text',color='black')
rc('axes',labelcolor='black',linewidth=2.0)
rc('font',size=15)
#rcParams["patch.force_edgecolor"] = True

import Calculate_ABA_Finkelsteinlab_Diewertje as CalcABA
#imp.reload(CalcABA);
import plotting_ABA_Diewertje as pltABA
#imp.reload(pltABA);
import plotting_Boyle_Diewertje as plt_B
#imp.reload(plt_B);

import sys 
sys.path.append('../code_general/')
#import CRISPR_free_energy_landscape as FreeEnergy
#imp.reload(FreeEnergy);
import read_model_ID;
# imp.reload(read_model_ID);

import sys 
sys.path.append('../code_general_Finkelsteinlab/')
import plotting_Finkelsteinlab as plt_F
#imp.reload(plt_F)


#import analysis_SA_fits as SAfits
#imp.reload(SAfits);


simset = []
no_good = []

chi_squared = [] 

#---------- collect simulations ---------------
for sim in range(1,13):
    sa = pd.read_csv('../Diewertje/1_3_2019/fit_1_3_2019_sim_' +str(sim)+'.txt', delimiter='\t', index_col=46)
    filename = '../Diewertje/1_3_2019/fit_1_3_2019_sim_' +str(sim) +'.txt'   
    chi_squared.append(sa.Potential.iloc[-1])
    simset.append(filename)
        
    
best_fit = simset[np.argmin(chi_squared)]
best_fit
# WRONG! SINCE SUCCESS IS FALSE!

# Load Parameters
model_id = 'init_limit_general_energies_v2'

# To find last line with the fitted parameters
f=open(best_fit)
lines=f.read().splitlines()
last_line=lines[-1]
last_line = last_line.split()
Finkel_params = list(map(float,last_line[:-2]))
print(Finkel_params)

# Nparams = (np.size(Finkel_params) #44
# boyle_params = plt_B.load_simm_anneal(filename, Nparams)

# Load data
IlyaData = pd.read_csv('../Data_ABA_Finkelsteinlab/cas9-target-e-replicate-1-delta-abas_Canonical_OT-r_0-2.csv')

plt.figure()
_ = pltABA.predict_single_mm(Finkel_params,model_id, data_file=IlyaData)
import numpy as np
from scipy import linalg
import sys
PATH_HPC05 = '/home/mklein1/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)
sys.path.append('../code_general/')
from read_model_ID import unpack_parameters
import Calculate_ABA_Finkelsteinlab as ABA


''''
Calculate the Chi-square function for Finkelstein data

Diewertje Dekker    Depken lab 
'''

def calc_Chi_square(parameters,xdata,ydata,yerr,
                    guide_length=20, model_id='general_energies'):
    # concentrations in nMolair
    # reference = 1 nMolair
    #mABA=ABA.calc_ABA(parameters, concentrations, reference, mismatch_positions, model_id = 'general_energies', guide_length = 20, T=10*60):
    Chi_square=0
    for i in range(len(xdata)):
        Chi_square=Chi_square+sum(((mABA[i]-ydata)/yerr)**2)
    return Chi_square
    


import numpy as np
from scipy import linalg
import sys
import Calculate_ABA_Finkelsteinlab_Diewertje as ABA
PATH_HPC05 = '/home/mklein1/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)
sys.path.append('../code_general/')
import read_model_ID


'''' 
Calculate the Chi-square function for Finkelstein data

Diewertje Dekker    Depken lab 
'''

def calc_Chi_square(parameters,xdata,ydata,yerr,concentrations,reference,ontarget_ABA,
                    guide_length=20,model_id='general_energies'):
    # concentrations in nMolair
    # reference = 1 nMolair
    mABA=ABA.calc_delta_ABA(parameters, concentrations, reference, xdata, ontarget_ABA, model_id =model_id, guide_length = guide_length, T=10*60)
    Chi_square=sum((((mABA*np.ones(len(ydata)))-np.array(ydata))/np.array(yerr))**2)
    return Chi_square
    


import numpy as np
from scipy import linalg
import sys
PATH_HPC05 = '/home/mklein1/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)
sys.path.append('../code_general/')
import read_model_ID
sys.path.append('../code_ABA_Finkelsteinlab/')
import Calculate_ABA_Finkelsteinlab as ABA


'''' 
Calculate the Chi-square function for Finkelstein data

Diewertje Dekker    Depken lab 
'''

def calc_Chi_square(parameters,xdata,ydata,yerr,concentrations,reference,
                    guide_length=20,model_id='general_energies'):
    # concentrations in nMolair
    # reference = 1 nMolair
    ontarget_ABA=ABA.calc_ABA(parameters,concentrations,reference, mismatch_positions=[],model_id = 'general_energies', guide_length = 20, T=10*60)
    Chi_square=0
    mABA=np.zeros(len(xdata))
    for i in range(len(xdata)):
        mABA[i]=ABA.calc_delta_ABA(parameters, concentrations, reference, xdata[i], ontarget_ABA, model_id = 'general_energies', guide_length = 20, T=10*60)
        Chi_square=Chi_square+sum(((mABA[i]-np.array(ydata[i]))/np.array(yerr[i]))**2)
    return Chi_square#, mABA
    


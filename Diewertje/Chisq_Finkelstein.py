import sys
# PATH_HPC05 = '/home/dddekker/BEP'
PATH_HPC05 = '/home/mklein1/Diewertje'
sys.path.append(PATH_HPC05)

import numpy as np
import Calculate_ABA_Finkelsteinlab_Diewertje as ABA


'''' 
Calculate the Chi-square function for Finkelstein data

Diewertje Dekker    Depken lab 
'''
# THis is for ABA in dataset! If you have delta ABA in dataset, change to: mABA=ABA.calc_ABA(...)

def calc_Chi_square(parameters,xdata,ydata,yerr,concentrations,reference,ontarget_ABA,
                    guide_length=20,model_id='general_energies'):
    # concentrations in nMolair
    # reference = 1 nMolair
    mABA=ABA.calc_ABA(parameters, concentrations, reference, xdata, model_id =model_id, guide_length = guide_length, T=10*60)
    Chi_square=sum((((mABA*np.ones(len(ydata)))-np.array(ydata))/np.array(yerr))**2)
    return Chi_square
    


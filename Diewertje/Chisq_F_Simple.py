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

def calc_Chi_square(parameters,xdata,ydata,yerr,concentrations,reference,ontarget_ABA,
                    guide_length=20,model_id='general_energies'):
    # concentrations in nMolair
    # reference = 1 nMolair
    epsilon = parameters[:-2]
    energies=ABA.get_energies(epsilon,xdata, guide_length)
    mABA=-np.log(np.sum(np.exp(-np.cumsum(energies))))
    Chi_square=sum((((mABA*np.ones(len(ydata)))-np.array(ydata))/np.array(yerr))**2)
    return Chi_square
    


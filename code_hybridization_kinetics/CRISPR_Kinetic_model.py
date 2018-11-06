import  numpy as np
import pandas as pd
import sys
sys.path.append('../code_general/')
sys.path.append('../code_Boyle/')

import read_model_ID
reload(read_model_ID)

from CRISPR_dCas9_binding_curve_Boyle import get_energies



'''
Implement Pclv from 'Hybridization Kinetics Explains CRISPR-Cas Off-Targeting Rules'. 

Allows us to use parameter sets from either binding or cleaving 
'''



def load_simm_anneal(filename, Nparams, fatch_solution='final'):
    '''
    Load the parameter set from simmulated annealing.
    Fix indexing based on parameterisation to correctly import the table
    :param filename:
    :param Nparams:
    :return:
    '''
    fit = pd.read_csv(filename, delimiter='\t', index_col=Nparams+2)
    fit = fit.reset_index()
    final_result = []
    for param in range(1, Nparams + 1):
        col = 'Parameter ' + str(param)

        if fatch_solution == 'final':
            final_result.append(fit[col].iloc[-1])
        else:
            final_result.append(fit[col].iloc[fatch_solution])

    sa_result = np.array(final_result)
    return sa_result

def get_backward_rates(energies, forward_rates,guide_length=20):
    '''
    Apply detailed balance condition to get the backward rates from the energies and forward rates

    :param energies:
    :param forward_rates:
    :param guide_length:
    :return:
    '''
    # 0) Construct array containing backward rates
    backward_rates = np.zeros(guide_length+1)

    # 1) Apply detailed balance condition:
    backward_rates = forward_rates[:-1] * np.exp(energies)
    return backward_rates


def translate_binding_to_cleaving(parameters, model_ID, rate_to_cleave, mismatch_positions,guide_length=20):
    '''
    Allows to use energy landscape inferred from binding data and converts it to parameter set used for cleavage.

    Binding parameters:
    1) Energies
    2) Forward rates

    Convert as follows:
    1) Use detailed balance to determine all the backward rates in the system:
        kb(n) = kf(n-1) * exp(+epsilon(n))
    2) Set an additional parameter: A constant catalytic rate from the complete R-loop
    (Automatically this will vary DeltaCLV is the final basepair is a mismatch)
    3) Use definition of kinetic biases to determine position of transition states:
        Delta(n) = - log(kb(n) /kf(n))

    Then you are ready to use these to deterine the transition landscape and calculate Pclv.
    :return:
    '''
    # 1) Use result from fit to binding data:
    epsilon, forward_rates = read_model_ID.unpack_parameters(parameters,model_ID,guide_length)

    energies = get_energies(epsilon,mismatch_positions, guide_length)


    # 2) Use detailed balance to get the backward rates:
    backward_rates = get_backward_rates(energies, forward_rates)

    # 3) Use Kramer's rate to convert to Delta's:
    forward_rates[-1] = rate_to_cleave
    # print 'backward: ', backward_rates[-1]
    # print 'forward: ', forward_rates[-1]
    Delta = np.log(forward_rates[1:]/backward_rates)
    return Delta


def get_transition_landscape(Delta):
    return np.cumsum(Delta)


def Pclv(Delta):
    DeltaT = get_transition_landscape(Delta)
    # For Pclv, I need the sum of the exponents of DeltaT_n
    exp_of_T = np.sum(np.exp(-DeltaT))
    return 1.0 / (1.0 + exp_of_T)




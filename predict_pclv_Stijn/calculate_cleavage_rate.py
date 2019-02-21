import numpy as np
import matplotlib.pylab as plt
from scipy import linalg
from scipy.optimize import curve_fit
import sys 
sys.path.append('../code_general/')
from read_model_ID import unpack_parameters



'''
Main functions
'''

def calc_chi_squared(parameters,mismatch_positions,ydata,yerr,times,
                    guide_length=20, model_id='Clv_init_limit_general_energies_v2'):
    
    k_model = calc_clv_rate(parameters, model_id, mismatch_positions, times,
                            guide_length)
    
    chi_sqrd = np.sum(((ydata-k_model)/yerr)**2)
    
    return chi_sqrd


def calc_clv_rate(parameters, model_id, mismatch_positions, times, guide_length=20):
    '''
    Calculates the cleavage rate, given the model_id, model parameters, 
    guide length and mismatch positions, at the given times.
    
    :param parameters:
    :param model_id:
    :param mismatch_positions:
    :param times:
    :param guide_length:
    :return: cleavage_rate
    '''
    #calculate master equation
    mat = get_master_equation(parameters=parameters, 
                              mismatch_positions=mismatch_positions, 
                              model_id=model_id, 
                              guide_length=guide_length)
    
    #calculate probabilities in time
    initial_prob = np.zeros(guide_length+2)
    initial_prob[0] = 1
    initial_prob.T
    
    prob_uncleaved = np.zeros(len(times))
    for i in range(len(times)):
        matrix_exponent = linalg.expm(+mat*times[i])
        prob_temp = matrix_exponent.dot(initial_prob)
        prob_uncleaved[i] = np.sum(prob_temp)
    
    
    #take the logarithm of the probabilities (zero values should not be considered)
    end = len(times)
           
    for i in range(len(times)):
        if prob_uncleaved[i]==0:
            end = i
            break
    
    times = times[0:end]
    prob_uncleaved = prob_uncleaved[0:end]
    
    prob_log = np.log(prob_uncleaved)
    
    #fit the log of the probability to a linear function, 
    #yielding the cleavage rate
    k, error = curve_fit(f,times,prob_log)
    return k[0]

'''
Linear fit-function used by calc_clv_rate
'''
def f(x,k):
    return -k*x




'''
Helper functions 
'''

def get_master_equation(parameters, mismatch_positions, model_id, guide_length):
    '''
    Construct rate matrix from given parameter set
    :param parameters:
    :param mismatch_positions:
    :param guide_length:
    :return:
    '''
    epsilon, forward_rates = unpack_parameters(parameters, model_id, guide_length)
    energies = get_energies(epsilon,mismatch_positions, guide_length)
    backward_rates = get_backward_rates(energies, forward_rates,guide_length )
    rate_matrix = build_rate_matrix(forward_rates, backward_rates)
    return rate_matrix

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
    epsI = new_epsilon[(guide_length+1):]
    energies = -1*new_epsilon[0:(guide_length+1)] # convention: epsC>0 means downward slope
    energies[0] = new_epsilon[0]                 # convention: epsPAM>0 means upward slope
    if len(mismatch_positions)>0:
        energies[mismatch_positions.astype(int)] += epsI[(mismatch_positions.astype(int)-1)]
    return energies


def get_backward_rates(energies, forward_rates,guide_length=20):
    '''
    Apply detailed balance condition to get the backward rates from the energies and forward rates

    :param energies:
    :param forward_rates:
    :param guide_length:
    :return:
    '''
    # 0) Construct array containing backward rates
    backward_rates = np.zeros(guide_length+2)

    # 1) Apply detailed balance condition:
    backward_rates[1:] = forward_rates[:-1] * np.exp(energies)

    # 2) No rate backward from solution state
    backward_rates[0] = 0.0
    return backward_rates


def build_rate_matrix(forward_rates, backward_rates):
    '''
    build matrix in Master Equation

    :param forward_rates:
    :param backward_rates:
    :return:
    '''
    diagonal1 = -(forward_rates + backward_rates)
    diagonal2 = backward_rates[1:]
    diagonal3 = forward_rates[:-1]
    # rate_matrix = np.zeros((len(forward_rates), len(forward_rates)))  # Build the matrix

    rate_matrix = np.diag(diagonal1, k=0) + np.diag(diagonal2, k=1) + np.diag(diagonal3, k=-1)

    # for diag in range(3):
    #     if diag == 0:
    #         rate_matrix = rate_matrix + np.diag(diagonal1, k=0)
    #     elif diag == 1:
    #         rate_matrix = rate_matrix + np.diag(diagonal2, k=1)
    #     elif diag == 2:
    #         rate_matrix = rate_matrix + np.diag(diagonal3, k=-1)

    return rate_matrix
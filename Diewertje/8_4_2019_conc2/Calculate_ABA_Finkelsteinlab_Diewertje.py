import numpy as np
from scipy import linalg
from scipy.optimize import curve_fit

import sys
# PATH_HPC05 = '/home/dddekker/BEP'
PATH_HPC05 = '/home/mklein1/Diewertje'
sys.path.append(PATH_HPC05)
sys.path.append('../code_general/') # added this after we run it on the cluster
from read_model_ID import unpack_parameters

def calc_Pbound(parameters, concentrations, reference, mismatch_positions, model_id = 'general_energies', guide_length = 20, T=10*60):
    #print('at the beginning', concentrations)
    rate_matrix = get_master_equation(parameters, mismatch_positions, model_id, guide_length)
    rel_concentration = concentrations/reference
    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))
    Pbound = []
    for c in rel_concentration:
        new_rate_matrix = rate_matrix.copy()
        new_rate_matrix[0][0] *= c
        new_rate_matrix[1][0] *= c
        Probability = get_Probability(new_rate_matrix, everything_unbound, T)
        Pbound.append(np.sum(Probability[1:]))
    Pbound = np.array(Pbound)
    concentrations = np.array(concentrations)
    conc3=True
    if(conc3==False):
        Kd, _ = curve_fit(Hill_eq, concentrations,Pbound,maxfev=10000)
    #print('at the end', concentrations)
    if(conc3==True):
        try: # because you do not always get a result from the fit
            Kd, _ = curve_fit(Hill_eq, concentrations,Pbound,maxfev=1000)
        except:
            Kd[0]= 0.0001 # 'error'
        if(Kd[0]<0): 
            Kd[0]=0.0001
    return Kd[0], Pbound, concentrations


def calc_ABA(parameters, concentrations, reference, mismatch_positions, model_id = 'general_energies', guide_length = 20, T=10*60):
    Kd, _,_ = calc_Pbound(parameters, concentrations, reference, mismatch_positions, model_id, guide_length, T)
    return np.log(Kd)

def calc_delta_ABA(parameters, concentrations, reference, mismatch_positions,ontarget_ABA, model_id = 'general_energies', guide_length = 20, T=10*60):
    ABA = calc_ABA(parameters, concentrations, reference, mismatch_positions, model_id, guide_length, T)
    return ABA-ontarget_ABA

def calc_Kd(parameters, concentrations, reference, mismatch_positions, model_id = 'general_energies', guide_length = 20, T=10*60):
    Kd, _,_ = calc_Pbound(parameters, concentrations, reference, mismatch_positions, model_id, guide_length, T)
    return Kd

'''
Helper functions 
'''

def Hill_eq(C, Kd):
    return (1.0+Kd/C)**(-1)

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
    epsI = np.array(new_epsilon[(guide_length+1):])
    energies = -1*np.array(new_epsilon[0:(guide_length+1)]) # convention: epsC>0 means downward slope
    energies[0] = new_epsilon[0]                 # convention: epsPAM>0 means upward slope
    if len(mismatch_positions)>0:
        energies[mismatch_positions.astype(int)] += epsI[(mismatch_positions.astype(int)-1)]
    energies=np.ndarray.tolist(energies)
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

    return rate_matrix

def get_Probability(rate_matrix, initial_condition,T=12*3600):
    '''
    solves the Master Equation for a given initial condition and desired time point
    :param rate_matrix: matrix with rates that makes up the Master Equation
    :param initial_condition: vector with initial configuration
    :param T: Evaluate solution at time T
    :return:
    '''
    P0 = initial_condition
    M = rate_matrix
    matrix_exponent = linalg.expm(+M*T)
    return matrix_exponent.dot(P0)

'''
Other useful functions 
'''
def calc_binding_curve(parameters, concentration, reference, mismatch_positions, T_list, model_id='general_energies', guide_length=20):

    rate_matrix = get_master_equation(parameters, mismatch_positions, model_id, guide_length)
    rel_concentration = concentration/reference
    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))
    new_rate_matrix = rate_matrix.copy()
    new_rate_matrix[0][0] *= rel_concentration
    new_rate_matrix[1][0] *= rel_concentration
    Probability = map(lambda t: 1-get_Probability(new_rate_matrix, everything_unbound, t)[0], T_list)
    return Probability

import numpy as np
import sys
sys.path.append('/home/svandersmagt/Energy_Landscape_dCas9/code_general/')
sys.path.append('../code_general/')
from scipy.optimize import curve_fit
from scipy import linalg
import read_model_ID as model
import calculate_cleavage_rate as clv


def coarse_grain_4state(parameters_full,model_id_full,mismatch=[],range_I=[7,13]):
    energies_full,rates_full = model.unpack_parameters(parameters_full,model_id_full)

    ## Mismatch
    for i in mismatch:
        energies_full[i] -= energies_full[i+20]

    ## energies20
    landscape = -np.cumsum(np.append(-energies_full[0],energies_full[1:21]))
    PAM = landscape[0] #energy PAM state
    Istate = np.argmin(landscape[range_I[0]:range_I[1]]) + range_I[0] #position of stable I state, within range_I
    energysum = 0.
    for i in range(range_I[0],range_I[1]):
        energysum += np.exp(-landscape[i])
    I = - np.log(energysum) #energy I state
    R = landscape[-1] #energy R state

    eps1 = PAM - I
    eps2 = I - R

    ## Getting matrix for barrier 1
    rates1 = rates_full[1:(Istate+1)]
    energies1 = -1*energies_full[1:(Istate+1)]
    backward1 = clv.get_backward_rates(energies1,rates1,Istate-2)
    matrix1 = clv.build_rate_matrix(rates1,backward1)

    ## Getting MFPT for barrier 1
    matrix1 *= -1
    Minv = np.linalg.inv(matrix1)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(Istate-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k1 = 1/MFPT

    ## Getting matrix for barrier 2
    rates2 = rates_full[(Istate+1):-1]
    energies2 = -1*energies_full[(Istate+1):-1]
    backward2 = clv.get_backward_rates(energies2,rates2,20-Istate-2)
    matrix2 = clv.build_rate_matrix(rates2,backward2)

    ## Getting MFPT for barrier 2
    matrix2 *= -1
    Minv = np.linalg.inv(matrix2)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(20-Istate-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k2 = 1/MFPT

    kb1 = np.exp(-eps1)*k1
    kb2 = np.exp(-eps2)*k2
    
    return [k1,k2,kb1,kb2],[-PAM,eps1,eps2],Istate

def calc_clv_rate_4state(rates,kon,kclv,PAM):
    
    forward_rates = np.array([kon,rates[0],rates[1],kclv])
    backward_rates = np.array([0.,kon*np.exp(-PAM),rates[2],rates[3]])
    matrix_coarse = clv.build_rate_matrix(forward_rates,backward_rates)
    matrix_coarse *= -1
    Minv = np.linalg.inv(matrix_coarse)
    vec = np.ones(len(Minv))
    everything_unbound = np.array([1.0] + [0.0]*(3))
    MFPT = vec.dot(Minv.dot(everything_unbound))
    k_coarse = 1/MFPT
    return k_coarse

def calc_aba_4state(rates,kon,PAM,concentrations,reference,T=10*60):
    
    forward_rates = np.array([kon,rates[0],rates[1],0.])
    backward_rates = np.array([0.,kon*np.exp(-PAM),rates[2],rates[3]])
    matrix_coarse = clv.build_rate_matrix(forward_rates,backward_rates)
    
    rel_concentration = concentrations/reference
    everything_unbound = np.array([1.0] + [0.0]*(3))
    Pbound = []
    for c in rel_concentration:
        new_rate_matrix = matrix_coarse.copy()
        new_rate_matrix[0][0] *= c
        new_rate_matrix[1][0] *= c
        Probability = get_Probability(new_rate_matrix, everything_unbound, T)
        Pbound.append(np.sum(Probability[1:]))
    Pbound = np.array(Pbound)
    concentrations = np.array(concentrations)
    conc3=True
    if(conc3==True):
        Kd, _ = curve_fit(Hill_eq, concentrations,Pbound,maxfev=1000)
        if(Kd[0]<0): 
            Kd[0]=0.0001
    return np.log(Kd[0])

def Hill_eq(C, Kd):
    return (1.0+Kd/C)**(-1)

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
    
def coarse_grain_4state_flexR(parameters_full,model_id_full,mismatch=[],range_I=[7,13],range_R=[16,21]):
    energies_full,rates_full = model.unpack_parameters(parameters_full,model_id_full)

    ## Mismatch
    for i in mismatch:
        energies_full[i] -= energies_full[i+20]

    ## energies20
    landscape = -np.cumsum(np.append(-energies_full[0],energies_full[1:21]))
    PAM = landscape[0] #energy PAM state
    Istate = np.argmin(landscape[range_I[0]:range_I[1]]) + range_I[0] #position of stable I state, within range_I
    energysum = 0.
    for i in range(range_I[0],range_I[1]):
        energysum += np.exp(-landscape[i])
    I = - np.log(energysum) #energy I state
    Rstate = np.argmin(landscape[range_R[0]:range_R[1]]) + range_R[0] #position of stable I state, within range_I
    energysum = 0.
    for i in range(range_R[0],range_R[1]):
        energysum += np.exp(-landscape[i])
    R = - np.log(energysum) #energy I state

    eps1 = PAM - I
    eps2 = I - R

    ## Getting matrix for barrier 1
    rates1 = rates_full[1:(Istate+1)]
    energies1 = -1*energies_full[1:(Istate+1)]
    backward1 = clv.get_backward_rates(energies1,rates1,Istate-2)
    matrix1 = clv.build_rate_matrix(rates1,backward1)

    ## Getting MFPT for barrier 1
    matrix1 *= -1
    Minv = np.linalg.inv(matrix1)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(Istate-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k1 = 1/MFPT
    
    ## Getting matrix for barrier 2
    rates2 = rates_full[(Istate+1):(Rstate+1)]
    energies2 = -1*energies_full[(Istate+1):(Rstate+1)]
    backward2 = clv.get_backward_rates(energies2,rates2,Rstate-Istate-2)
    matrix2 = clv.build_rate_matrix(rates2,backward2)

    ## Getting MFPT for barrier 2
    matrix2 *= -1
    Minv = np.linalg.inv(matrix2)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(Rstate-Istate-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k2 = 1/MFPT

    kb1 = np.exp(-eps1)*k1
    kb2 = np.exp(-eps2)*k2
    
    kclv_effective = np.exp(-landscape[-1])/np.exp(-R)*rates_full[-1]
    
    return [k1,k2,kb1,kb2],[-PAM,eps1,eps2],Istate,Rstate,kclv_effective

def coarse_grain_4state_coarseClv(parameters_full,model_id_full,mismatch=[],range_I=[7,13],range_R=[16,21]):
    energies_full,rates_full = model.unpack_parameters(parameters_full,model_id_full)

    ## Mismatch
    for i in mismatch:
        energies_full[i] -= energies_full[i+20]

    ## energies20
    landscape = -np.cumsum(np.append(-energies_full[0],energies_full[1:21]))
    PAM = landscape[0] #energy PAM state
    Istate = np.argmin(landscape[range_I[0]:range_I[1]]) + range_I[0] #position of stable I state, within range_I
    energysum = 0.
    for i in range(range_I[0],range_I[1]):
        energysum += np.exp(-landscape[i])
    I = - np.log(energysum) #energy I state
    Rstate = np.argmin(landscape[range_R[0]:range_R[1]]) + range_R[0] #position of stable I state, within range_I
    energysum = 0.
    for i in range(range_R[0],range_R[1]):
        energysum += np.exp(-landscape[i])
    R = - np.log(energysum) #energy I state

    eps1 = PAM - I
    eps2 = I - R

    ## Getting matrix for barrier 1
    rates1 = rates_full[1:(Istate+1)]
    energies1 = -1*energies_full[1:(Istate+1)]
    backward1 = clv.get_backward_rates(energies1,rates1,Istate-2)
    matrix1 = clv.build_rate_matrix(rates1,backward1)

    ## Getting MFPT for barrier 1
    matrix1 *= -1
    Minv = np.linalg.inv(matrix1)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(Istate-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k1 = 1/MFPT
    
    ## Getting matrix for barrier 2
    rates2 = rates_full[(Istate+1):(Rstate+1)]
    energies2 = -1*energies_full[(Istate+1):(Rstate+1)]
    backward2 = clv.get_backward_rates(energies2,rates2,Rstate-Istate-2)
    matrix2 = clv.build_rate_matrix(rates2,backward2)

    ## Getting MFPT for barrier 2
    matrix2 *= -1
    Minv = np.linalg.inv(matrix2)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(Rstate-Istate-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k2 = 1/MFPT

    kb1 = np.exp(-eps1)*k1
    kb2 = np.exp(-eps2)*k2
    
    ## Getting matrix for kclv
    rates3 = rates_full[(Rstate+1):]
    energies3 = -1*energies_full[(Rstate+1):]
    backward3 = clv.get_backward_rates(energies3,rates3,20-Rstate-1)
    matrix3 = clv.build_rate_matrix(rates3,backward3)

    ## Getting MFPT for kclv
    matrix3 *= -1
    Minv = np.linalg.inv(matrix3)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(20-Rstate))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    kclv_effective = 1/MFPT
    
    return [k1,k2,kb1,kb2],[-PAM,eps1,eps2],Istate,Rstate,kclv_effective
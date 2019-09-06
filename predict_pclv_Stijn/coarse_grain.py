import numpy as np
import sys
sys.path.append('/home/svandersmagt/Energy_Landscape_dCas9/code_general/')
sys.path.append('../code_general/')
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
    everything_atPAM = np.array([1.0] + [0.0]*(3))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k_coarse = 1/MFPT
    return k_coarse

def coarse_grain_5state(parameters_full,model_id_full,mismatch=[],range_I1=[2,7],range_I2=[7,13]):
    energies_full,rates_full = model.unpack_parameters(parameters_full,model_id_full)

    ## Mismatch
    for i in mismatch:
        energies_full[i] -= energies_full[i+20]

    ## energies20
    landscape = -np.cumsum(np.append(-energies_full[0],energies_full[1:21]))
    PAM = landscape[0] #energy PAM state
    Istate1 = np.argmin(landscape[range_I1[0]:range_I1[1]]) + range_I1[0] #position of stable I1 state, within range_I1
    Istate2 = np.argmin(landscape[range_I2[0]:range_I2[1]]) + range_I2[0] #position of stable I2 state, within range_I2
    energysum = 0.
    for i in range(range_I1[0],range_I1[1]):
        energysum += np.exp(-landscape[i])
    I1 = - np.log(energysum) #energy I1 state
    energysum = 0.
    for i in range(range_I2[0],range_I2[1]):
        energysum += np.exp(-landscape[i])
    I2 = - np.log(energysum) #energy I2 state
    R = landscape[-1] #energy R state

    eps1 = PAM - I1
    eps2 = I1 - I2
    eps3 = I2 - R

    ## Getting matrix for barrier 1
    rates1 = rates_full[1:(Istate1+1)]
    energies1 = -1*energies_full[1:(Istate1+1)]
    backward1 = clv.get_backward_rates(energies1,rates1,Istate1-2)
    matrix1 = clv.build_rate_matrix(rates1,backward1)

    ## Getting MFPT for barrier 1
    matrix1 *= -1
    Minv = np.linalg.inv(matrix1)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(Istate1-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k1 = 1/MFPT

    ## Getting matrix for barrier 2
    rates2 = rates_full[(Istate1+1):(Istate2+1)]
    energies2 = -1*energies_full[(Istate1+1):(Istate2+1)]
    backward2 = clv.get_backward_rates(energies2,rates2,(Istate2-Istate1)-2)
    matrix2 = clv.build_rate_matrix(rates2,backward2)

    ## Getting MFPT for barrier 2
    matrix2 *= -1
    Minv = np.linalg.inv(matrix2)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(Istate2-Istate1-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k2 = 1/MFPT
    
    ## Getting matrix for barrier 3
    rates3 = rates_full[(Istate2+1):-1]
    energies3 = -1*energies_full[(Istate2+1):-1]
    backward3 = clv.get_backward_rates(energies3,rates3,20-Istate2-2)
    matrix3 = clv.build_rate_matrix(rates3,backward3)

    ## Getting MFPT for barrier 3
    matrix3 *= -1
    Minv = np.linalg.inv(matrix3)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(20-Istate2-1))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k3 = 1/MFPT

    kb1 = np.exp(-eps1)*k1
    kb2 = np.exp(-eps2)*k2
    kb3 = np.exp(-eps3)*k3
    
    return [k1,k2,k3,kb1,kb2,kb3],[-PAM,eps1,eps2,eps3],[Istate1,Istate2]

def calc_clv_rate_5state(rates,kon,kclv,PAM):
    
    forward_rates = np.array([kon,rates[0],rates[1],rates[2],kclv])
    backward_rates = np.array([0.,kon*np.exp(-PAM),rates[3],rates[4],rates[5]])
    matrix_coarse = clv.build_rate_matrix(forward_rates,backward_rates)
    matrix_coarse *= -1
    Minv = np.linalg.inv(matrix_coarse)
    vec = np.ones(len(Minv))
    everything_atPAM = np.array([1.0] + [0.0]*(4))
    MFPT = vec.dot(Minv.dot(everything_atPAM))
    k_coarse = 1/MFPT
    return k_coarse
    
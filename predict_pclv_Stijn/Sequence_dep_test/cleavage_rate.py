import numpy as np
from scipy import linalg
from scipy.optimize import curve_fit
import sys 
sys.path.append('/home/svandersmagt/Energy_Landscape_dCas9/code_general/')
sys.path.append('../../code_general/')
from read_model_ID import unpack_parameters



'''
Main functions
'''

def calc_chi_squared(parameters,xdata,ydata,yerr,
                    guide_length, model_id):
       
    k_model = calc_clv_rate_fast(parameters, model_id, xdata,
                            guide_length)
        
    ydata = np.array(ydata)
    yerr = np.array(yerr)
    chi_sqrd = np.sum(((ydata-np.log10(k_model))/yerr)**2)
    return chi_sqrd


def calc_clv_rate(parameters, model_id, xdata, guide_length=20):
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
     times = np.array([0.0,12.0,60.0,180.0,600.0,1800.0,6000.0,18000.0,60000.0])
     
     '''
     Linear fit-function used by calc_clv_rate
     '''
     def f(x,k):
         return -k*x
     
     
     #calculate master equation
     mat = get_master_equation(parameters=parameters, 
                       xdata=xdata, 
                       model_id=model_id, 
                       guide_length=guide_length)
     
     #calculate probabilities in time
     initial_prob = np.zeros(guide_length+2)
     initial_prob[0] = 1
     
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



def calc_clv_rate_fast(parameters, model_id, xdata, guide_length=20):
    #calculate master equation
    mat = get_master_equation(parameters=parameters, 
                              xdata=xdata, 
                              model_id=model_id, 
                              guide_length=guide_length)
    
    M = -1 * mat
    try:
        Minv = np.linalg.inv(M)
    except:
        print 'INVERTING MATRIX FAILED, USE SLOWER METHOD'
        sys.stdout.flush()
        return calc_clv_rate(parameters, model_id, xdata)
    vec = np.ones(len(Minv))
    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))
    MFPT = vec.dot(Minv.dot(everything_unbound))
    k = 1/MFPT
   
    return k

    
'''
Helper functions 
'''

def get_master_equation(parameters, xdata, model_id, guide_length):
    '''
    Construct rate matrix from given parameter set
    :param parameters:
    :param mismatch_positions:
    :param guide_length:
    :return:
    '''
    epsilonConfig, epsilonBind, forward_rates = unpack_parameters(parameters, model_id, guide_length)
    energies = get_energies(epsilonConfig,epsilonBind,xdata, guide_length)
    backward_rates = get_backward_rates(energies, forward_rates,guide_length )
    rate_matrix = build_rate_matrix(forward_rates, backward_rates)
    return rate_matrix

def get_energies(epsilonConfig,epsilonBind,xdata, guide_length=20):
    '''
    For general (position dependent) model make a list with the energies at every bound state
    At positions with a mismatch incorporated: add mismatch penalty (epsI[mm_pos])

    So this function returns the minima in the energy lanscape (the actual energy at every state)

    :param epsilon: [epsPAM, epsC[state_1],....,epsC[state_N],epsI[state_1],...epsI[state_N] ]
    provide as np.array()
    :param mismatch_positions: each mismach position has a range [1,2, ... , 20]
    :return: vector with the minima in the landscape
    '''
    
    energies = -1*epsilonConfig # convention: epsC>0 means downward slope
    energies[0] = epsilonConfig[0] # convention: epsPAM>0 means upward slope
    guide = xdata[1]
    target = xdata[0]
    for pos in range(guide_length):
        if guide[pos]=='A':
            if target[pos]=='T':
                energies[pos+1]+=epsilonBind[0]
            elif target[pos]=='A':
                energies[pos+1]+=epsilonBind[1]
            elif target[pos]=='C':
                energies[pos+1]+=epsilonBind[2]
            elif target[pos]=='G':
                energies[pos+1]+=epsilonBind[3]
        elif guide[pos]=='T':
            if target[pos]=='T':
                energies[pos+1]+=epsilonBind[4]
            elif target[pos]=='A':
                energies[pos+1]+=epsilonBind[7]
            elif target[pos]=='C':
                energies[pos+1]+=epsilonBind[10]
            elif target[pos]=='G':
                energies[pos+1]+=epsilonBind[12]
        elif guide[pos]=='G':
            if target[pos]=='T':
                energies[pos+1]+=epsilonBind[2]
            elif target[pos]=='A':
                energies[pos+1]+=epsilonBind[5]
            elif target[pos]=='C':
                energies[pos+1]+=epsilonBind[8]
            elif target[pos]=='G':
                energies[pos+1]+=epsilonBind[9]
        elif guide[pos]=='C':
            if target[pos]=='T':
                energies[pos+1]+=epsilonBind[3]
            elif target[pos]=='A':
                energies[pos+1]+=epsilonBind[6]
            elif target[pos]=='C':
                energies[pos+1]+=epsilonBind[9]
            elif target[pos]=='G':
                energies[pos+1]+=epsilonBind[11]
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

def calc_association_rate(rate_matrix,timepoints=[500.,1000.,1500.],guide_length=20,rel_concentration=0.1):
    '''
    Association experiment for the effective on-rate:
    (repeating the protocol as used in Boyle et al.)
    :return:
    '''

    #1) Calculate bound fraction for specified time points:
    # Association experiment starts with all dCas9 being unbound:
    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))

    #3) Association rate is taken at 1nM the other data at 10nM --> adjust k_on appropriatly:
    new_rate_matrix = rate_matrix.copy()
    new_rate_matrix[0][0] *= rel_concentration  #rel_concentration=1 corresponds to 10nM
    new_rate_matrix[1][0] *= rel_concentration

    bound_fraction = []
    for time in timepoints:
        Probabilities = get_Probability(rate_matrix=new_rate_matrix,initial_condition=everything_unbound,T=time)
        bound_fraction.append( np.sum(Probabilities[1:]))

    #2) Use least-squares to fit straight line with origin forced through zero:
    kon_lin_fit = least_squares_line_through_origin(x_points=np.array(timepoints),y_points=np.array(bound_fraction))
    return kon_lin_fit

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

def least_squares_line_through_origin(x_points, y_points):
    return np.sum( x_points*y_points )/np.sum( x_points*x_points )
import numpy as np
from scipy import linalg
from scipy.optimize import curve_fit
import sys 
sys.path.append('/home/svandersmagt/Energy_Landscape_dCas9/code_general/')
sys.path.append('../code_general/')
import read_model_ID as model


'''
Main functions
'''

def calc_chi_squared(parameters,mismatch_positions,ydata,yerr,chi_weights,combined_fit,
                    guide_length, model_id):
    
    if not combined_fit:
        k_model = calc_clv_rate_fast(parameters, model_id, mismatch_positions,
                                guide_length)
            
        ydata = np.array(ydata)
        yerr = np.array(yerr)
        
        chi_sqrd_clv_perfect = 0.0
        chi_sqrd_clv_single = 0.0
        chi_sqrd_clv_double = 0.0
        
        if len(mismatch_positions)==0:
            chi_sqrd_clv_perfect = np.sum(((ydata-np.log10(k_model))/yerr)**2)
        
        elif len(mismatch_positions)==1:
            chi_sqrd_clv_single = np.sum(((ydata-np.log10(k_model))/yerr)**2)
        
        elif len(mismatch_positions)==2:
            chi_sqrd_clv_double = np.sum(((ydata-np.log10(k_model))/yerr)**2)
            
        chi_sqrd = (chi_sqrd_clv_perfect*chi_weights[0] +
                    chi_sqrd_clv_single*chi_weights[1] +
                    chi_sqrd_clv_double*chi_weights[2])
        
        return chi_sqrd

    if combined_fit:
        k_model_clv, k_model_on = calc_clv_on(parameters, model_id,
                                              mismatch_positions, guide_length)
        ydata_clv = np.array(ydata[0])
        ydata_on = np.array(ydata[1])
        yerr_clv = np.array(yerr[0])
        yerr_on = np.array(yerr[1])
        
        chi_sqrd_clv_perfect = 0.0
        chi_sqrd_clv_single = 0.0
        chi_sqrd_clv_double = 0.0
        chi_sqrd_on_perfect = 0.0
        chi_sqrd_on_single = 0.0
        chi_sqrd_on_double = 0.0
        
        if len(mismatch_positions)==0:
            chi_sqrd_clv_perfect = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
            chi_sqrd_on_perfect = np.sum(((ydata_on-k_model_on)/yerr_on)**2)
        
        elif len(mismatch_positions)==1:
            chi_sqrd_clv_single = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
            chi_sqrd_on_single = np.sum(((ydata_on-k_model_on)/yerr_on)**2)
        
        elif len(mismatch_positions)==2:
            chi_sqrd_clv_double = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
            chi_sqrd_on_double = np.sum(((ydata_on-k_model_on)/yerr_on)**2)
            
        chi_sqrd = (chi_sqrd_clv_perfect*chi_weights[0] + chi_sqrd_on_perfect*chi_weights[3] +
                    chi_sqrd_clv_single*chi_weights[1] + chi_sqrd_on_single*chi_weights[4] +
                    chi_sqrd_clv_double*chi_weights[2] + chi_sqrd_on_double*chi_weights[5])
        
        return chi_sqrd

def calc_clv_rate(parameters, model_id, mismatch_positions, guide_length=20):
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
                       mismatch_positions=mismatch_positions, 
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



def calc_clv_rate_fast(parameters, model_id, mismatch_positions, guide_length=20):
    #calculate master equation
    mat = get_master_equation(parameters=parameters, 
                              mismatch_positions=mismatch_positions, 
                              model_id=model_id, 
                              guide_length=guide_length)
    
    M = -1 * mat
    try:
        Minv = np.linalg.inv(M)
    except:
        print 'INVERTING MATRIX FAILED, USE SLOWER METHOD'
        sys.stdout.flush()
        return calc_clv_rate(parameters, model_id, mismatch_positions)
    vec = np.ones(len(Minv))
    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))
    MFPT = vec.dot(Minv.dot(everything_unbound))
    k = 1/MFPT
   
    return k


def calc_clv_on(parameters, model_id, mismatch_positions, guide_length):
    
    matrix_clv, matrix_on = get_master_equation_clv_on(parameters,
                                                       mismatch_positions,
                                                       model_id,
                                                       guide_length)
    
    ## Calculating cleavage rate
    M_clv = -1 * matrix_clv
    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))
    try:
        Minv_clv = np.linalg.inv(M_clv)
        vec_clv = np.ones(len(Minv_clv))
        MFPT_clv = vec_clv.dot(Minv_clv.dot(everything_unbound))
        k_clv = 1/MFPT_clv
    except:
        print 'INVERTING MATRIX FAILED, USE SLOWER METHOD'
        sys.stdout.flush()
        if len(parameters)==44:
            parameters_clv = parameters[1:42] + parameters[42:44]
        if len(parameters)==43:
            parameters_clv = parameters[0:40] + parameters[41:43]
        k_clv = calc_clv_rate(parameters_clv, model_id[0], mismatch_positions, guide_length)
    
    
    ## Calculating on rate
    k_on = calc_association_rate(rate_matrix=matrix_on,timepoints=[500.,1000.,1500.],
                                                 guide_length=guide_length)
    
    return k_clv, k_on
    
'''
Helper functions 
'''

def get_master_equation_clv_on(parameters,mismatch_positions,model_id,guide_length):
    
    model_id_clv,model_id_on,parameters_clv,parameters_on = model.combined_model(parameters,model_id)
    
    epsilon_on, forward_rates_on = model.unpack_parameters(parameters_on, model_id_on, guide_length)
    energies_on = get_energies(epsilon_on,mismatch_positions, guide_length)
    backward_rates_on = get_backward_rates(energies_on, forward_rates_on,guide_length )
    matrix_on = build_rate_matrix(forward_rates_on, backward_rates_on)
    
    epsilon_clv, forward_rates_clv = model.unpack_parameters(parameters_clv, model_id_clv, guide_length)
    energies_clv = get_energies(epsilon_clv,mismatch_positions, guide_length)
    backward_rates_clv = get_backward_rates(energies_clv,forward_rates_clv,guide_length)
    
    #most of the matrix entries are equal
    matrix_clv = matrix_on.copy()
    
    #different rate from solution
    matrix_clv[0][0] = -forward_rates_clv[0]
    matrix_clv[1][0] = forward_rates_clv[0]
    
    #different cleavage rate
    matrix_clv[-1][-1] = -(forward_rates_clv[-1]+backward_rates_clv[-1]) 
    
    #different ePAM
    matrix_clv[0][1] = backward_rates_clv[1]
    matrix_clv[1][1] = -(backward_rates_clv[1]+forward_rates_clv[1])
    matrix_clv[1][2] = backward_rates_clv[2]
    matrix_clv[2][2] = -(backward_rates_clv[2]+forward_rates_clv[2])
    
    return matrix_clv, matrix_on

def get_master_equation(parameters, mismatch_positions, model_id, guide_length):
    '''
    Construct rate matrix from given parameter set
    :param parameters:
    :param mismatch_positions:
    :param guide_length:
    :return:
    '''
    epsilon, forward_rates = model.unpack_parameters(parameters, model_id, guide_length)
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
    #backward_rates[1:] = forward_rates[:-1] * np.exp(energies), somehow this suddenly does not work anymore..
    for i in range(1,len(backward_rates)):
        backward_rates[i] = forward_rates[i-1]*np.exp(energies[i-1])
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
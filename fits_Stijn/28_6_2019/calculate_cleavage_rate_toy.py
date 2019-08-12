import numpy as np
from scipy import linalg
from scipy.optimize import curve_fit
import read_model_ID_toy as model


'''
Main functions
'''

def calc_chi_squared(parameters,mismatch_positions,ydata,yerr,chi_weights,combined_fit,
                     model_id,log_on=False):
    
    if not combined_fit:
        k_model = calc_cleavage_rate_fast(parameters, model_id, mismatch_positions)
        if k_model < 10**-5:
            k_model = 10**-5
            
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
                                              mismatch_positions)
        if k_model_clv < 10**-5:
            k_model_clv = 10**-5
        
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
        
        
        
        if not log_on:
            if len(mismatch_positions)==0:
                chi_sqrd_clv_perfect = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
                chi_sqrd_on_perfect = np.sum(((ydata_on-k_model_on)/yerr_on)**2)

            elif len(mismatch_positions)==1:
                chi_sqrd_clv_single = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
                chi_sqrd_on_single = np.sum(((ydata_on-k_model_on)/yerr_on)**2)

            elif len(mismatch_positions)==2:
                chi_sqrd_clv_double = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
                chi_sqrd_on_double = np.sum(((ydata_on-k_model_on)/yerr_on)**2)
        
        else:
            if len(mismatch_positions)==0:
                chi_sqrd_clv_perfect = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
                chi_sqrd_on_perfect = np.sum(((ydata_on-np.log10(k_model_on))/yerr_on)**2)

            elif len(mismatch_positions)==1:
                chi_sqrd_clv_single = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
                chi_sqrd_on_single = np.sum(((ydata_on-np.log10(k_model_on))/yerr_on)**2)

            elif len(mismatch_positions)==2:
                chi_sqrd_clv_double = np.sum(((ydata_clv-np.log10(k_model_clv))/yerr_clv)**2)
                chi_sqrd_on_double = np.sum(((ydata_on-np.log10(k_model_on))/yerr_on)**2)
            
        chi_sqrd = (chi_sqrd_clv_perfect*chi_weights[0] + chi_sqrd_on_perfect*chi_weights[3] +
                    chi_sqrd_clv_single*chi_weights[1] + chi_sqrd_on_single*chi_weights[4] +
                    chi_sqrd_clv_double*chi_weights[2] + chi_sqrd_on_double*chi_weights[5])
        
        return chi_sqrd

def calc_cleavage_rate_fast(parameters,model_id,mismatch_positions):
    
    rate_matrix = master_equation(parameters,model_id,mismatch_positions)
    
    M = -1 * rate_matrix
    try:
        Minv = np.linalg.inv(M)
    except:
        print 'INVERTING MATRIX FAILED, USE SLOWER METHOD'
        sys.stdout.flush()
        return calc_clv_rate_slow(parameters, model_id, mismatch_positions)
    vec = np.ones(len(Minv))
    everything_unbound = np.array([1.0] + [0.0] * (len(rate_matrix[0])-1))
    MFPT = vec.dot(Minv.dot(everything_unbound))
    k = 1/MFPT
    
    return k

def calc_cleavage_rate_slow(parameters,model_id,mismatch_positions):
    
    times = np.array([0.0,12.0,60.0,180.0,600.0,1800.0,6000.0,18000.0,60000.0])
     
    def f(x,k):
        return -k*x
     
     
    #calculate master equation
    mat = master_equation(parameters,model_id,mismatch_positions)
     
    #calculate probabilities in time
    initial_prob = np.zeros(len(mat[0]))
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
    
def master_equation(parameters,model_id,mismatch_positions):
    forward_combined, backward_combined = combined_rates(parameters,model_id,mismatch_positions)
    
    diagonal1 = -(forward_combined + backward_combined)
    diagonal2 = backward_combined[1:]
    diagonal3 = forward_combined[:-1]
    #rate_matrix = np.zeros((len(forward_rates), len(forward_rates)))  # Build the matrix

    rate_matrix = np.diag(diagonal1, k=0) + np.diag(diagonal2, k=1) + np.diag(diagonal3, k=-1)
    #rate_matrix[-2][-2] -= forward_combined[-1]/100. # to add cleavage from the second to last state
    
    return rate_matrix

def combined_rates(parameters,model_id,mismatch_positions):
    epsilon,rates = model.unpack_parameters(parameters,model_id)
    if type(mismatch_positions)==type([]):
            mismatch_positions = np.array(mismatch_positions)
    
    #energies
    new_epsilon = epsilon.copy()        
    energies_OT = -1*new_epsilon[0:3]
    energies_OT[0] *= -1 #PAM, I, 0, R, 0
    barrier1 = new_epsilon[3]
    barrier2 = new_epsilon[4]
    
    #On target rates
    rates = rates.copy() #rate from solution and cleavage rate remain
    rates[1] *= np.exp(-barrier1) #rate over first barrier
    rates[2] *= np.exp(-barrier2) #rate over second barrier
    
    backward = np.zeros(len(rates))
    backward[0] = 0.0 #no rate backward from solution
    backward[1] = rates[0]*np.exp(energies_OT[0]) #backward rate from PAM
    backward[2] = rates[1]*np.exp(energies_OT[1]) #backward rate over first barrier
    backward[3] = rates[2]*np.exp(energies_OT[2]) #backward rate over second barrier

    #handling mismatches
    for i in mismatch_positions:
        if i < 8: #mismatch in the first barrier, could be placed at 8 or 9 or 10. difference can be compensated by mismatch penalty 8,9
            rates[1] *= np.exp(-new_epsilon[4+i])
        elif i < 12: #mismatch in first flat part, pretty sure it should be at 12
            backward[2] *= np.exp(new_epsilon[4+i])
        elif i < 18: #mismatch in second barrier, does not matter that much I think(maybe it does for the engineered proteins)
            rates[2] *= np.exp(-new_epsilon[4+i])
        elif i <21: #mismatch in last flat part
            backward[3] *= np.exp(new_epsilon[4+i])

    return rates, backward

def calc_clv_on(parameters, model_id, mismatch_positions):
    
    matrix_clv, matrix_on = get_master_equation_clv_on(parameters,
                                                       mismatch_positions,
                                                       model_id)
    
    ## Calculating cleavage rate
    M_clv = -1 * matrix_clv
    everything_unbound = np.array([1.0] + [0.0] * (len(matrix_clv[0])-1))
#    try:
    Minv_clv = np.linalg.inv(M_clv)
    vec_clv = np.ones(len(Minv_clv))
    MFPT_clv = vec_clv.dot(Minv_clv.dot(everything_unbound))
    k_clv = 1/MFPT_clv
#    except:
#        print 'INVERTING MATRIX FAILED, USE SLOWER METHOD'
#        sys.stdout.flush()
#        if len(parameters)==44:
#            parameters_clv = parameters[1:42] + parameters[42:44]
#        if len(parameters)==43:
#            parameters_clv = parameters[0:40] + parameters[41:43]
#        k_clv = calc_clv_rate(parameters_clv, model_id[0], mismatch_positions, guide_length)
    
    
    ## Calculating on rate
    k_on = calc_association_rate(rate_matrix=matrix_on,timepoints=[500.,1000.,1500.])
    
    return k_clv, k_on

def get_master_equation_clv_on(parameters,mismatch_positions,model_id):
    
    model_id_clv,model_id_on,parameters_clv,parameters_on = model.combined_model(parameters,model_id)
    
    matrix_on = master_equation(parameters_on,model_id_on,mismatch_positions)
    
    forward_rates_clv,backward_rates_clv = combined_rates(parameters_clv,model_id_clv,mismatch_positions)
    
    #most of the matrix entries are equal
    matrix_clv = matrix_on.copy()
    
    #different rate from solution
    matrix_clv[0][0] = -forward_rates_clv[0]
    matrix_clv[1][0] = forward_rates_clv[0]
    
    #different cleavage rate
    matrix_clv[-1][-1] = -(forward_rates_clv[-1]+backward_rates_clv[-1])
    #matrix_clv[-2][-2] -= forward_rates_clv[-1]/100. # to add cleavage from the second to last state
    
    #different ePAM
    matrix_clv[0][1] = backward_rates_clv[1]
    matrix_clv[1][1] = -(backward_rates_clv[1]+forward_rates_clv[1])
    matrix_clv[1][2] = backward_rates_clv[2]
    matrix_clv[2][2] = -(backward_rates_clv[2]+forward_rates_clv[2])
    
    return matrix_clv, matrix_on


def calc_association_rate(rate_matrix,timepoints=[500.,1000.,1500.],rel_concentration=0.1):
    '''
    Association experiment for the effective on-rate:
    (repeating the protocol as used in Boyle et al.)
    :return:
    '''

    #1) Calculate bound fraction for specified time points:
    # Association experiment starts with all dCas9 being unbound:
    everything_unbound = np.array([1.0] + [0.0] * (len(rate_matrix[0])-1))

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
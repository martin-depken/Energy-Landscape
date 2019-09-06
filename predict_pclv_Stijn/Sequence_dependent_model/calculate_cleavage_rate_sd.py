import numpy as np
from scipy import linalg
from scipy.optimize import curve_fit
import sys 
sys.path.append('/home/svandersmagt/Energy_Landscape_dCas9/code_general/')
sys.path.append('../../code_general/')
from read_model_ID_sd import unpack_parameters



'''
Main functions
'''

def calc_chi_squared(parameters,xdata,ydata,yerr,
                    guide_length, model_id,before,after):
       
    k_model = calc_clv_rate_fast(parameters, model_id, xdata,
                            guide_length,before,after)
        
    ydata = np.array(ydata)
    yerr = np.array(yerr)
    chi_sqrd = np.sum(((ydata-np.log10(k_model))/yerr)**2)
    return chi_sqrd


def calc_clv_rate(parameters, model_id, xdata, guide_length=20, before=True, after=False):
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
                       guide_length=guide_length, before=before, after=after)
     
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



def calc_clv_rate_fast(parameters, model_id, xdata, guide_length=20, before=True, after=False):
    #calculate master equation
    mat = get_master_equation(parameters=parameters, 
                              xdata=xdata, 
                              model_id=model_id, 
                              guide_length=guide_length, before=before, after=after)
    
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

def get_master_equation(parameters, xdata, model_id, guide_length, before, after):
    '''
    Construct rate matrix from given parameter set
    :param parameters:
    :param mismatch_positions:
    :param guide_length:
    :return:
    '''
    epsilon_si,rates_si,epsilon_i_sd = unpack_parameters(parameters, model_id, guide_length)
    energies = get_energies(epsilon_si,epsilon_i_sd,xdata, guide_length, before, after)
    backward_rates = get_backward_rates(energies, rates_si,guide_length )
    rate_matrix = build_rate_matrix(rates_si, backward_rates)
    return rate_matrix

def get_energies(epsilon_si,epsilon_i_sd,xdata, guide_length, before, after):
    '''
    For general (position dependent) model make a list with the energies at every bound state
    At positions with a mismatch incorporated: add mismatch penalty (epsI[mm_pos])

    So this function returns the minima in the energy lanscape (the actual energy at every state)

    :param epsilon: [epsPAM, epsC[state_1],....,epsC[state_N],epsI[state_1],...epsI[state_N] ]
    provide as np.array()
    :param mismatch_positions: each mismach position has a range [1,2, ... , 20]
    :return: vector with the minima in the landscape
    '''
    
    ## Sequence independent part:
    mismatch_positions = xdata[2]
    if type(mismatch_positions)==type([]):
        mismatch_positions = np.array(mismatch_positions)
    
    new_epsilon = epsilon_si.copy()
    epsilon_i_si = new_epsilon[(guide_length+1):]
    energies = -1*new_epsilon[0:(guide_length+1)] # convention: epsC>0 means downward slope
    energies[0] = new_epsilon[0]                 # convention: epsPAM>0 means upward slope
    if len(mismatch_positions)>0:
        energies[mismatch_positions.astype(int)] += epsilon_i_si[(mismatch_positions.astype(int)-1)]
    
    if before and not after: #neighbour before position counts
        ## Sequence dependent part:
        if len(mismatch_positions)>0:
            guide = xdata[1]
            target = xdata[0]
            for pos in mismatch_positions:
                pos -= 1
                if pos == 0:
                    ## guide[pos] == C at this position
                    if target[pos]=='T':
                        energies[pos+1]+=epsilon_i_sd[26] #P-CA
                    if target[pos]=='A':
                        energies[pos+1]+=epsilon_i_sd[29] #P-CT
                    if target[pos]=='G':
                        energies[pos+1]+=epsilon_i_sd[32] #P-CC
                else:
                    if guide[pos]=='A':
                        if target[pos]=='T':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[0] #A-AA
                            elif guide[pos-1]=='T':
                                energies[pos+1]+=epsilon_i_sd[1] #U-AA
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[2] #G-AA
                            elif guide[pos-1]=='C':
                                energies[pos-1]+=epsilon_i_sd[3] #C-AA

                        elif target[pos]=='C':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[4] #A-AG
                            elif guide[pos-1]=='T':
                                energies[pos+1]+=epsilon_i_sd[5] #U-AG
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[6] #G-AG
                            elif guide[pos-1]=='C':
                                energies[pos-1]+=epsilon_i_sd[7] #C-AG

                        elif target[pos]=='G':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[8] #A-AC
                            elif guide[pos-1]=='T':
                                energies[pos+1]+=epsilon_i_sd[9] #U-AC
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[10] #G-AC
                            elif guide[pos-1]=='C':
                                energies[pos-1]+=epsilon_i_sd[11] #C-AC



                    elif guide[pos]=='T':
                        if target[pos]=='A':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[12] #A-UT
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[13] #G-UT

                        if target[pos]=='C':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[14] #A-UG
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[15] #G-UG

                        if target[pos]=='G':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[16] #A-UC
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[17] #G-UC



                    elif guide[pos]=='G':
                        if target[pos]=='T':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[18] #A-GA
                            elif guide[pos-1]=='C':
                                energies[pos+1]+=epsilon_i_sd[19] #C-GA

                        if target[pos]=='A':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[20] #A-GT
                            elif guide[pos-1]=='C':
                                energies[pos+1]+=epsilon_i_sd[21] #C-GT

                        if target[pos]=='C':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[22] #A-GG
                            elif guide[pos-1]=='C':
                                energies[pos+1]+=epsilon_i_sd[23] #C-GG



                    elif guide[pos]=='C':
                        if target[pos]=='T':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[24] #A-CA
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[25] #G-CA

                        if target[pos]=='A':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[27] #A-CT
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[28] #G-CT

                        if target[pos]=='G':
                            if guide[pos-1]=='A':
                                energies[pos+1]+=epsilon_i_sd[30] #A-CC
                            elif guide[pos-1]=='G':
                                energies[pos+1]+=epsilon_i_sd[31] #G-CC
                                
    if after and not before: #neighbour after position counts:
        if len(mismatch_positions)>0:
            guide = xdata[1]
            target = xdata[0]
            for pos in mismatch_positions:
                pos -= 1
                if pos == 19:
                    ## guide[pos] == G at this position
                    if target[pos]=='T':
                        energies[pos+1]+=epsilon_i_sd[15] #GA-A
                    if target[pos]=='A':
                        energies[pos+1]+=epsilon_i_sd[18] #GT-A
                    if target[pos]=='C':
                        energies[pos+1]+=epsilon_i_sd[21] #GG-A
                else:
                    if guide[pos]=='A':
                        if target[pos]=='T':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[0] #AA-A
                            elif guide[pos+1]=='T':
                                energies[pos+1]+=epsilon_i_sd[1] #AA-U
                            elif guide[pos+1]=='G':
                                energies[pos+1]+=epsilon_i_sd[2] #AA-G
                            elif guide[pos+1]=='C':
                                energies[pos-1]+=epsilon_i_sd[3] #AA-C

                        elif target[pos]=='C':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[4] #AG-A
                            elif guide[pos+1]=='T':
                                energies[pos+1]+=epsilon_i_sd[5] #AG-U
                            elif guide[pos+1]=='G':
                                energies[pos+1]+=epsilon_i_sd[6] #AG-G
                            elif guide[pos+1]=='C':
                                energies[pos-1]+=epsilon_i_sd[7] #AG-C

                        elif target[pos]=='G':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[8] #AC-A
                            elif guide[pos+1]=='T':
                                energies[pos+1]+=epsilon_i_sd[9] #AC-U
                            elif guide[pos+1]=='G':
                                energies[pos+1]+=epsilon_i_sd[10] #AC-G
                            elif guide[pos+1]=='C':
                                energies[pos-1]+=epsilon_i_sd[11] #AC-C



                    elif guide[pos]=='T':
                        if target[pos]=='A':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[12] #UT-A

                        if target[pos]=='C':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[13] #UG-A


                        if target[pos]=='G':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[14] #UC-A




                    elif guide[pos]=='G':
                        if target[pos]=='T':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[15] #GA-A
                            elif guide[pos+1]=='T':
                                energies[pos+1]+=epsilon_i_sd[16] #GA-U
                            elif guide[pos+1]=='C':
                                energies[pos+1]+=epsilon_i_sd[17] #GA-C

                        if target[pos]=='A':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[18] #GT-A
                            elif guide[pos+1]=='T':
                                energies[pos+1]+=epsilon_i_sd[19] #GT-U
                            elif guide[pos+1]=='C':
                                energies[pos+1]+=epsilon_i_sd[20] #GT-C

                        if target[pos]=='C':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[21] #GG-A
                            elif guide[pos+1]=='T':
                                energies[pos+1]+=epsilon_i_sd[22] #GG-U
                            elif guide[pos+1]=='C':
                                energies[pos+1]+=epsilon_i_sd[23] #GG-C



                    elif guide[pos]=='C':
                        if target[pos]=='T':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[24] #CA-A
                            elif guide[pos+1]=='G':
                                energies[pos+1]+=epsilon_i_sd[25] #CA-G

                        if target[pos]=='A':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[26] #CT-A
                            elif guide[pos+1]=='G':
                                energies[pos+1]+=epsilon_i_sd[27] #CT-G

                        if target[pos]=='G':
                            if guide[pos+1]=='A':
                                energies[pos+1]+=epsilon_i_sd[28] #CC-A
                            elif guide[pos+1]=='G':
                                energies[pos+1]+=epsilon_i_sd[29] #CC-G
                                
                      
    if before and after: #both neighbours
        ## Sequence dependent part:
        if len(mismatch_positions)>0:
            guide = xdata[1]
            target = xdata[0]
            for pos in mismatch_positions:
                pos -= 1
                if pos == 0:
                    ## guide[pos] == P-C-G at this position
                    if target[pos]=='T':
                        energies[pos+1]+=epsilon_i_sd[38] #P-CA-G
                    if target[pos]=='A':
                        energies[pos+1]+=epsilon_i_sd[41] #P-CT-G
                    if target[pos]=='G':
                        energies[pos+1]+=epsilon_i_sd[44] #P-CC-G
                elif pos == 19:
                    ## guide[pos] == A-G-A at this position
                    if target[pos]=='T':
                        energies[pos+1]+=epsilon_i_sd[27] #A-GA-A
                    if target[pos]=='A':
                        energies[pos+1]+=epsilon_i_sd[30] #A-GT-A
                    if target[pos]=='C':
                        energies[pos+1]+=epsilon_i_sd[33] #A-GG-A
                else:
                    if guide[pos-1]=='C' and guide[pos]=='G' and guide[pos+1]=='C':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[29] #C-GA-C
                        elif target[pos]=='A':
                            energies[pos+1]+=epsilon_i_sd[32] #C-GT-C
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[35] #C-GG-C
                            
                    elif guide[pos-1]=='G' and guide[pos]=='C' and guide[pos+1]=='A':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[37] #G-CA-A
                        elif target[pos]=='A':
                            energies[pos+1]+=epsilon_i_sd[40] #G-CT-A
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[43] #G-CC-A
                            
                    elif guide[pos-1]=='C' and guide[pos]=='A' and guide[pos+1]=='G':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[6] #C-AA-G
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[13] #C-AG-G
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[20] #C-AC-G
                            
                    elif guide[pos-1]=='A' and guide[pos]=='G' and guide[pos+1]=='A':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[27] #A-GA-A
                        elif target[pos]=='A':
                            energies[pos+1]+=epsilon_i_sd[30] #A-GT-A
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[33] #A-GG-A
                            
                    elif guide[pos-1]=='G' and guide[pos]=='A' and guide[pos+1]=='G':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[5] #G-AA-G
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[12] #G-AG-G
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[19] #G-AC-G
                            
                    elif guide[pos-1]=='A' and guide[pos]=='G' and guide[pos+1]=='T':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[28] #A-GA-U
                        elif target[pos]=='A':
                            energies[pos+1]+=epsilon_i_sd[31] #A-GT-U
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[34] #A-GG-U
                            
                    elif guide[pos-1]=='G' and guide[pos]=='T' and guide[pos+1]=='A':
                        if target[pos]=='A':
                            energies[pos+1]+=epsilon_i_sd[22] #G-UT-A
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[24] #G-UG-A
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[26] #G-UC-A
                            
                    elif guide[pos-1]=='T' and guide[pos]=='A' and guide[pos+1]=='G':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[2] #U-AA-G
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[9] #U-AG-G
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[16] #U-AC-G
                            
                    elif guide[pos-1]=='G' and guide[pos]=='A' and guide[pos+1]=='A':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[4] #G-AA-A
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[11] #G-AG-A
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[18] #G-AC-A
                            
                    elif guide[pos-1]=='A' and guide[pos]=='A' and guide[pos+1]=='A':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[0] #A-AA-A
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[7] #A-AG-A
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[14] #A-AC-A
                            
                    elif guide[pos-1]=='A' and guide[pos]=='A' and guide[pos+1]=='T':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[1] #A-AA-U
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[8] #A-AG-U
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[15] #A-AC-U
                            
                    elif guide[pos-1]=='A' and guide[pos]=='T' and guide[pos+1]=='A':
                        if target[pos]=='A':
                            energies[pos+1]+=epsilon_i_sd[21] #A-UT-A
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[23] #A-UG-A
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[25] #A-UC-A
                            
                    elif guide[pos-1]=='T' and guide[pos]=='A' and guide[pos+1]=='C':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[3] #U-AA-C
                        elif target[pos]=='C':
                            energies[pos+1]+=epsilon_i_sd[10] #U-AG-C
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[17] #U-AC-C
                            
                    elif guide[pos-1]=='A' and guide[pos]=='C' and guide[pos+1]=='G':
                        if target[pos]=='T':
                            energies[pos+1]+=epsilon_i_sd[36] #A-CA-G
                        elif target[pos]=='A':
                            energies[pos+1]+=epsilon_i_sd[39] #A-CT-G
                        elif target[pos]=='G':
                            energies[pos+1]+=epsilon_i_sd[42] #A-CC-G
                            
                    else:
                        print 'huh wat dom'


        

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

def calc_ABA(parameters, concentrations, reference, xdata, model_id = 'general_energies', guide_length = 20, T=10*60, before=True, after=False):
    Kd, _,_ = calc_Pbound(parameters, concentrations, reference, xdata, model_id, guide_length, T, before, after)
    return np.log(Kd)

def calc_Pbound(parameters, concentrations, reference, xdata, model_id = 'general_energies', guide_length = 20, T=10*60, before=True, after=False):
    #print('at the beginning', concentrations)
    rate_matrix = get_master_equation(parameters, xdata, model_id, guide_length,before,after)
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
    conc3=False
    if(conc3==False):
        Kd, _ = curve_fit(Hill_eq, concentrations,Pbound,maxfev=10000)
    return Kd[0], Pbound, concentrations

def Hill_eq(C, Kd):
    return (1.0+Kd/C)**(-1)
'''
Code to compare multiple fits from Simmulated Annealing together.

Can get the distributions/boxplots of different parameters amongst the different runs

Also built-in a sweep to find the bounds for the internal forward rate based on the averaged (free-)energy landscape
'''
import importlib as imp
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('ticks');
current_colors = sns.color_palette()

import sys
sys.path.append('../code_general/')
from read_model_ID import unpack_parameters
import CRISPR_free_energy_landscape_Diewertje as FreeEnergy
imp.reload(FreeEnergy);
sys.path.append('../code_Boyle/')
import plotting_Boyle_Diewertje as plt_B
import CRISPR_dCas9_binding_curve_Boyle as dCas9
imp.reload(dCas9);

from scipy import optimize
import Boyle_data_processing as Bdata
imp.reload(Bdata);

#import CRISPR_dCas9_binding_curve_Boyle as dCas9
#reload(dCas9);
sys.path.append('../code_Pclv/')
import CRISPR_Kinetic_model_Diewertje as Pclv
imp.reload(Pclv);

import get_parameters_fit_Diewertje as gpf

########################

def P(kon,koff,t):
    return kon/(kon+koff) * (1 - np.exp(-(kon+koff)*t))

def least_squares_line_through_origin(x_points, y_points):
    return np.sum( x_points*y_points )/np.sum( x_points*x_points )


def association(kon,koff):
    times = np.array([0.5,1.0,1.5])
    fluoresence = P(kon,koff,times)
    k_asso = least_squares_line_through_origin(x_points=times, y_points=fluoresence)
    return k_asso

def find_on_target(kon,k_asso):
    return k_asso - association(kon,koff=0.0)

def calc_on_rate():
    '''
    using the measured association rate for the on-target, invert the relation between association and on-rate
    through root-finding.

    uses 'brentq()' algorithm from scipy to find the root(s).

    :return: estimated on-target on rate
    '''

    BoyleData = Bdata.read()
    OnTarget = BoyleData[BoyleData.MM_pos.apply(len)==0]
    asso_rate_OT = float(OnTarget.on_slope)*1000
    kOT = optimize.brentq(find_on_target,a=0.01, b=1.0,args=(asso_rate_OT))
    return kOT/1000.


def process_SA_fits(simset,Nparams=44, model_id='init_limit_general_energies_v2'):
    '''
    Load the desired fits and store into dataframes:
    1. the value of epsilon_C at every position
    2. the value of epsilon_I at every position
    3. the value of the energy_landscape at every position
    4. the bound free-energy for the on-target at every position
    5. the estimated on-rate for the on-target. Both with and without assuming the solution and PAM have equillibrated

    :param simset:     Let every entry from simset be a filename including path

    :return: 6 separate pandas dataframes
    '''


    # prepare some stuff for the data frames:
    column_names = ['sim', 'PAM']
    for n in range(1, 21):
        column_names.append(str(n))

    match_rows = []
    mismatch_rows = []
    rates_rows = []
    landscape_rows = []
    FreeEnergy_rows = []
    fast_Rloop_rows = []

    for filename in simset:
        parameters = gpf.load_simm_anneal(filename, Nparams)
        parameters[0]=3.5 # If we want to fix epsilon PAM
        epsilon, fwd_rates = unpack_parameters(parameters, model_id)
        Energy_landscape, FreeEnergy_landscape = FreeEnergy.plot_free_energy_landscape(parameters, model_id,
                                                                                       show_plot=False);

        match_data = {}
        match_data['sim'] = filename
        for i, key in enumerate(column_names[1:]):
            match_data[key] = epsilon[i]
        match_rows.append(match_data)

        mismatch_data = {}
        mismatch_data['sim'] = filename
        for i, key in enumerate(column_names[2:]):
            mismatch_data[key] = epsilon[21 + i]
        mismatch_rows.append(mismatch_data)

        rates_data = {}
        rates_data['sim'] = filename
        rates_data['R1_to_PAM'] = fwd_rates[1] * np.exp(-epsilon[1])
        rates_data['PAM_to_sol'] = fwd_rates[0] * np.exp(epsilon[0])
        for i, key in enumerate(['sol_to_PAM', 'PAM_to_R1', 'internal']):
            rates_data[key] = fwd_rates[i]
        rates_rows.append(rates_data)

        landscape_data = {}
        landscape_data['sim'] = filename
        for i, key in enumerate(['sol'] + column_names[1:]):
            landscape_data[key] = Energy_landscape[i]
        landscape_rows.append(landscape_data)

        FreeEnergy_data = {}
        FreeEnergy_data['sim'] = filename
        for i, key in enumerate(column_names[1:]):
            FreeEnergy_data[key] = FreeEnergy_landscape[i]
        FreeEnergy_rows.append(FreeEnergy_data)

        rate_to_cleave = 10 ** 9
        Delta = Pclv.translate_binding_to_cleaving(parameters, model_id, rate_to_cleave, mismatch_positions=[])
        P = Pclv.Pclv(Delta[1:])
        kPR = fwd_rates[1] * P
        kSP = 0.1 * fwd_rates[0]
        kPS = kSP * np.exp(epsilon[0] + np.log(10.0))
        k_OT = kPR * kSP / (kPR + kSP + kPS)
        kinetic = k_OT * 1000

        # Assume PAM and solution equillibrates: Initiating the R-loop is limmiting
        thermodynamicPAM = kPR / (1.0 + np.exp(+epsilon[0] + np.log(10.0))) * 1000

        # Assume that PAM and (first) R-loop state equillibrate: Binding from solution is limmiting


        P2 = Pclv.Pclv(Delta[2:])
        kf = fwd_rates[2]
        E1 = epsilon[1]
        k_OT =  (kSP * kf*P2/np.exp(E1))/(kSP + kPS + kf*P2/np.exp(E1))
        thermodynamicR  =  k_OT*1000

        fast_Rloop_data = {}
        fast_Rloop_data['sim'] = filename
        fast_Rloop_data['kinetic'] = kinetic
        fast_Rloop_data['eq_PAM'] = thermodynamicPAM
        fast_Rloop_data['eq_PR'] = thermodynamicR
        fast_Rloop_rows.append(fast_Rloop_data)

    matches = pd.DataFrame(match_rows, columns=column_names)
    mismatches = pd.DataFrame(mismatch_rows, columns=['sim'] + [str(i) for i in range(1, 21)])
    rates = pd.DataFrame(rates_rows, columns=['sim', 'sol_to_PAM', 'PAM_to_sol', 'PAM_to_R1', 'R1_to_PAM', 'internal'])
    landscape = pd.DataFrame(landscape_rows, columns=['sim', 'sol'] + column_names[1:])
    free_energy = pd.DataFrame(FreeEnergy_rows, columns=column_names)
    fast_Rloop = pd.DataFrame(fast_Rloop_rows, columns=['sim', 'kinetic', 'eq_PAM','eq_PR'])

    matches.set_index('sim', inplace=True)
    mismatches.set_index('sim', inplace=True)
    rates.set_index('sim', inplace=True)
    landscape.set_index('sim', inplace=True)
    free_energy.set_index('sim', inplace=True)
    fast_Rloop.set_index('sim', inplace=True)
    return matches, mismatches, rates, landscape, free_energy, fast_Rloop


def average_solution(simset):
    '''
    use some identifiers to choose the simulations I want to average over

    to easily toggle between sims in different folders maybe make a variable 'sim_set'
    and use:
    for sim in sim_set:
        load simulation
        ...
    '''
    matches, mismatches, rates, landscape, free_energy, fast_Rloop = process_SA_fits(simset)

    # 1) get average energy landscape:
    landscape_avg = np.array(landscape.mean())

    # 2) get average mismatch penalties:
    epsI_avg = np.array(mismatches.mean())

    # 3) get average on-target binding rate:
    kOT_avg = float(fast_Rloop.mean()['kinetic'])

    return landscape_avg, epsI_avg, kOT_avg/1000.0


def median_solution(simset):
    '''
    use some identifiers to choose the simulations I want to average over

    to easily toggle between sims in different folders maybe make a variable 'sim_set'
    and use:
    for sim in sim_set:
        load simulation
        ...
    '''
    matches, mismatches, rates, landscape, free_energy, fast_Rloop = process_SA_fits(simset)

    # 1) get average energy landscape:
    landscape_avg = np.array(landscape.median())

    # 2) get average mismatch penalties:
    epsI_avg = np.array(mismatches.median())

    # 3) get average on-target binding rate:
    kOT_avg = float(fast_Rloop.median()['kinetic'])
    # kOT_avg = float(fast_Rloop.median()['thermodynamic'])

    return landscape_avg, epsI_avg, kOT_avg/1000.0




def calc_rate_PAM_to_Rloop(energy_landscape, kf, kOT, kSP=1000.0):
    '''
    Get kPR from kOT and the two other fitted forward rates.

    Adjusted to use exact equation involving both kf and kSP
    :param energy_landscape:
    :param kf:
    :param kOT:
    :param kSP:
    :return:
    '''


    Epsilon = np.diff(energy_landscape)
    Delta = -np.diff(energy_landscape)
    Delta[0] *= -1

    E_SP = Epsilon[0]
    # ---- associaton rate data is taken at 1nM, we report the parameter values at 10nM ----
    kSP *= 0.1
    Kd = np.exp(+E_SP + np.log(10.0))

    P2 = Pclv.Pclv(Delta[2:])
    c = kOT * (1 + Kd)  /(1 - kOT/kSP)
    E_RP = Epsilon[1]
    alpha = 1.0 /(kf * P2) * np.exp(+E_RP) * c
    kPR = c * (1 - alpha) ** (-1)
    return kPR , alpha


def test_rates(kf, kOT, energy_landscape, epsilon_I,
               xdata, ydata, yerr, model_id='init_limit_general_energies_v2',
               kSP=1000.0):
    '''
    calculate the chi-squared value based on the single-mismatch association rate data
    for a set of particular values of the internal forward rate and the rate from solution to PAM
    '''

    # 1) Use kOT, the energy landscape, kf and kSP to determine kPR:
    kPR, _ = calc_rate_PAM_to_Rloop(energy_landscape, kf, kOT, kSP)

    # 2) Calculate chi-squared with single mismatch off-targets:
    Epsilon = np.diff(energy_landscape)
    Epsilon[1:] *= -1

    new_parameters = list(Epsilon) + list(epsilon_I)

    new_parameters.append(np.log10(kSP))
    new_parameters.append(np.log10(kPR))
    new_parameters.append(np.log10(kf))
    new_parameters = np.array(new_parameters)

    ONtarget_occupancy = 1.0

    V = 0
    for i in range(len(xdata)):
        V += dCas9.calc_Chi_square(new_parameters, xdata[i], ydata[i], yerr[i], ONtarget_occupancy,
                                   model_id=model_id)

    return V, new_parameters


def optimize_internal_forward_rate(simset, forward_rates, mode='mean',kSP=1000.0):
    xdata, ydata, yerr = Bdata.prepare_multiprocessing(use_single_mm_only=True,
                                                       use_on_rate=True,
                                                       use_off_rate=False,
                                                       use_occupancy=False)

    if mode == 'mean':
        landscape_avg, epsI_avg, kOT = average_solution(simset)
    elif mode =='median':
        landscape_avg, epsI_avg, kOT = median_solution(simset)

    # kOT = calc_on_rate()
    _, lower_bnd_kf = calc_rate_PAM_to_Rloop(kSP=kSP, kf=1.0, kOT=kOT, energy_landscape=landscape_avg)
    V = []
    kf_vals = []
    for kf in forward_rates:
        if kf < lower_bnd_kf:
            continue
        kf_vals.append(kf)
        v, _ = test_rates(kf=kf, kSP=kSP,
                                 kOT=kOT,
                                 energy_landscape=landscape_avg,
                                 epsilon_I=epsI_avg,
                                 xdata=xdata,
                                 ydata=ydata,
                                 yerr=yerr)
        V.append(v)


    # Now find the optimum:
    kf_opt = kf_vals[np.argmin(V)]
    V_opt, parameters_opt = test_rates(kf=kf_opt,
                                       kSP=kSP,
                                      kOT=kOT,
                                      energy_landscape=landscape_avg,
                                      epsilon_I=epsI_avg,
                                      xdata=xdata,
                                      ydata=ydata,
                                      yerr=yerr)

    return V, kf_vals, kf_opt, parameters_opt, V_opt


def grid_search_forward_rates(simset, int_forward_rates, sol_to_PAM_rates,
                              save_to_file = True,
                              mode='median',
                              today='30/11/2018'):
    '''
    2D grid search to find both internal forward rate and rate from solution to PAM based on the average landscape
    '''

    # --- for every value of kSP, perform a 1D grid search along kf ----

    grid_V = np.nan * np.ones((len(int_forward_rates), len(sol_to_PAM_rates)))
    grid_kfvals = np.nan * np.ones((len(int_forward_rates), len(sol_to_PAM_rates)))

    # partial optimum along kf:
    grid_kfopt = []
    grid_Vopt = []
    grid_parameters = []
    for i, kSP in enumerate(sol_to_PAM_rates):
        V, kf_vals, kf_opt, parameters_opt, V_opt = optimize_internal_forward_rate(simset, int_forward_rates, mode,
                                                                                   kSP=kSP)

        grid_V[-len(V):, i] = V
        grid_kfvals[-len(kf_vals):, i] = kf_vals
        grid_kfopt.append(kf_opt)
        grid_Vopt.append(V_opt)
        grid_parameters.append(parameters_opt)

    # --- find optimum on grid -----
    Vopt = np.min(grid_Vopt)
    kf_opt = grid_kfopt[np.argmin(grid_Vopt)]
    kSP_opt = sol_to_PAM_rates[np.argmin(grid_Vopt)]
    parameters_opt = grid_parameters[np.argmin(grid_Vopt)]

    if save_to_file:
        #---- store parameters into text file ----
        fit_info = str(len(simset)) + ' in total ,  folder 25_10_2018/sims: 1 - 150 & 19_10_2018'
        model_id = 'init_limit_general_energies_v2'
        file_params = '../data/25_10_2018/' + mode + '_landscape_Boyle_2Dgrid.txt'
        write_parameters(parameters_opt, model_id, file_params, today, fit_info, mode)

        # ---- store grid points into Excel file ------
        filename = '../data/25_10_2018/grid_search_' + mode +  today.replace('/','_') +'.xlsx'
        write_grid(grid_V, grid_kfvals, Vopt, kf_opt, kSP_opt, int_forward_rates, sol_to_PAM_rates, filename)


    return grid_V, grid_kfvals, Vopt, kf_opt, kSP_opt, parameters_opt


def write_parameters(parameters, model_id, filename,
                    today, fit_info, mode):

    O = open(filename, 'w')
    O.write('# date of modification: ' + today + '\n')
    O.write('# model ID: '+ model_id + '\n')
    O.write('# SA fits used: ' + fit_info + '\n' )
    O.write('# analysis: '+ mode + '\n')

    for param in parameters:
        O.write(str(param) + '\n')
    O.close()
    return


def write_grid(V, kf_vals, Vopt, kf_opt, kSP_opt, forward_rates,binding_rates, filename):
    # ------ calculated chi-squared values on grid points ------------
    df = pd.DataFrame(V)
    convert_col_names = {}
    for i in range(len(binding_rates)):
        convert_col_names[i] = binding_rates[i]

    convert_row_names = {}
    for j in range(len(forward_rates)):
        convert_row_names[j] = forward_rates[j]
    df.rename(index=convert_row_names, columns=convert_col_names, inplace=True)

    # --- accepted forward rates -------
    df2 = pd.DataFrame(kf_vals)
    df2.rename(index=convert_row_names, columns=convert_col_names, inplace=True)

    # ----- optimal set of kf and kSP on the grid ------
    df3 = pd.DataFrame()
    df3['kf'] = [kf_opt]
    df3['kSP'] = [kSP_opt]
    df3['V'] = [Vopt]

    # ----- save it all into an excel file -----
    ExcelFile = pd.ExcelWriter(path=filename)
    df.to_excel(ExcelFile, sheet_name='chi_squared')
    df2.to_excel(ExcelFile, sheet_name='forward_rates')
    df3.to_excel(ExcelFile, sheet_name='optimum')
    ExcelFile.save()
    ExcelFile.close()
    return
















###############################################################
def replace_lower_triangle(M):
    '''
    sets the values in the lower triangle of matrix M equal to 1.
    The model predictions contain both mismatch config (i,j) and (j,i).
    Setting half of the association rates equal to 1.0 will prevent those duplicates
    from contributing to the selection done below
    :param M:
    :return:
    '''
    K = M
    for i in range(len(M)):
        for j in range(len(M)):
            if i > j:
                K[i, j] = 1
    return K



def difference_model_predictions(model, reference):
    '''
    Calculate average distance between two model predictions per datapoint
    :param model:
    :param reference:
    :return:
    '''
    reference = replace_lower_triangle(reference)
    model = replace_lower_triangle(model)
    N = len(reference)
    total_nmbr_of_points = N * (N + 1) * 0.5
    sum_difference = np.sum(np.abs(model - reference) / (reference))
    diff = sum_difference/total_nmbr_of_points
    return diff



def select_on_prediction(simset, chi_squared, percentage,
                         Nparams=44,
                         model_id='init_limit_general_energies_v2',
                         precalculated=False, score=None,
                         save_scores=True, filename='select_with_predcitions.txt'
                         ):
    '''
    Select those solutions that whose model prediction on the training data differs no more than x% from
    the prediction belonging to the solution with the lowest chi-squared value in the set of simulations.
    '''
    if not precalculated:
        # ----- Start selection: Retrieve the best fit first --------
        chi_squared = np.array(chi_squared)
        best_fit = simset[np.argmin(chi_squared)]
        parameters = gpf.load_simm_anneal(best_fit, Nparams)
        _, model_best, _ = plt_B.calc_predictions(parameters=parameters, model_id=model_id)
        # model_best = replace_lower_triangle(model_best)

        # ----- Compare difference in model prediction to this best fit ---
        score = []
        for sim in simset:
            parameters = gpf.load_simm_anneal(sim, Nparams)
            _, model, _ = plt_B.calc_predictions(parameters=parameters, model_id=model_id)
            diff = difference_model_predictions(model, model_best)
            score.append(  diff )
        score = np.array(score)

    # ----- select simulations whose difference in predicted values differs less then x% from the best fit ----
    selected_scores = score[score <= percentage]

    simset = np.array(simset)
    selected_sims = simset[score <= percentage]

    # ---- return selected_sims and selected_scores ------
    if save_scores:
        np.savetxt(filename, score)

    return selected_sims, selected_scores, score












def select_on_chi2(Chi2, simset, percentage=0.05):
    '''
    Select those solutions that whose chi-squared differs no more than x% from
    the lowest chi-squared value in the set of simulations.
    '''
    # ---- find minimum value of Chi-squared amongst replicates ---
    best_solution = min(Chi2)

    # ---- select those simulations with a chi2 that differs no more than x% from the best solution ---
    selected_solutions = Chi2[Chi2 <= ((1 + percentage) * best_solution)]

    # ----- select the simulations files ----
    simset = np.array(simset)
    selected_sims = simset[Chi2 <= ((1 + percentage) * best_solution)]
    return selected_solutions, selected_sims, (1 + percentage) * best_solution

def Tukey_outlier_test(data, simset, k=1.5):
    '''
    Use 'Tukey outlier test' to filter outliers from dataset

    protocol:
    1. Assume the data to be Gaussian distributed (null-hypothesis)
    2. A datapoint is NOT an outlier if it lies within the 25% and 75% percentiles
    3. p-value is the probability of obtaining a datapoint outside this interval. Prob. to reject the null-hypothsis.
    '''

    # ---- Tukey's test says: x in [Q1-1.5*IQR, Q3+1.5*IQR], with IQR=Q3-Q1 ----
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    low = Q1 - k * IQR
    high = Q3 + k * IQR

    # ---- Choose the datapoints ---
    Tukey_data = data[(data <= high) & (data >= low)]

    # ---- select corresponding simulations -----
    simset = np.array(simset)
    Tukey_sims = simset[(data <= high) & (data >= low)]
    return Tukey_data, Tukey_sims, low, high



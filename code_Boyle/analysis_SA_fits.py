'''
Code to compare multiple fits from Simmulated Annealing together.

Can get the distributions/boxplots of different parameters amongst the different runs

Also built-in a sweep to find the bounds for the internal forward rate based on the averaged (free-)energy landscape
'''

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('ticks');
current_colors = sns.color_palette()

import sys
sys.path.append('../code_general/')
from read_model_ID import unpack_parameters
import CRISPR_free_energy_landscape as FreeEnergy
reload(FreeEnergy);
import plotting_Boyle as plt_B
reload(plt_B)
import CRISPR_dCas9_binding_curve_Boyle as dCas9
reload(dCas9);

from scipy import optimize
import Boyle_data_processing as Bdata
reload(Bdata);

import CRISPR_dCas9_binding_curve_Boyle as dCas9
reload(dCas9);
sys.path.append('../code_Pclv/')
import CRISPR_Kinetic_model as Pclv
reload(Pclv);

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
        parameters = plt_B.load_simm_anneal(filename, Nparams)
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
        k_init = kPR * kSP / (kPR + kSP + kPS)
        kinetic = k_init * 1000
        thermodynamic = kPR / (1.0 + np.exp(+epsilon[0] + np.log(10.0))) * 1000

        fast_Rloop_data = {}
        fast_Rloop_data['sim'] = filename
        fast_Rloop_data['kinetic'] = kinetic
        fast_Rloop_data['thermodynamic'] = thermodynamic
        fast_Rloop_rows.append(fast_Rloop_data)

    matches = pd.DataFrame(match_rows, columns=column_names)
    mismatches = pd.DataFrame(mismatch_rows, columns=['sim'] + [str(i) for i in range(1, 21)])
    rates = pd.DataFrame(rates_rows, columns=['sim', 'sol_to_PAM', 'PAM_to_sol', 'PAM_to_R1', 'R1_to_PAM', 'internal'])
    landscape = pd.DataFrame(landscape_rows, columns=['sim', 'sol'] + column_names[1:])
    free_energy = pd.DataFrame(FreeEnergy_rows, columns=column_names)
    fast_Rloop = pd.DataFrame(fast_Rloop_rows, columns=['sim', 'kinetic', 'thermodynamic'])

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




def calc_rate_PAM_to_Rloop(kf, kOT, energy_landscape):
    Epsilon = np.diff(energy_landscape)
    Delta = -np.diff(energy_landscape)
    Delta[0] *= -1

    Kd = np.exp(+Epsilon[0] + np.log(10.0))

    P2 = Pclv.Pclv(Delta[2:])
    c = kOT * (1 + Kd)
    E1 = Epsilon[1]
    alpha = 1.0 /(kf * P2) * np.exp(+E1) * c
    kPR = c * (1 - alpha) ** (-1)
    return kPR , alpha

def test_forward_rate(kf, kOT, energy_landscape, epsilon_I,
                     xdata, ydata, yerr, model_id='init_limit_general_energies_v2'):
    '''
    calculate the chi-squared value based on the single-mismatch association rate data
    for a particular value of the internal forward rate
    '''

    # 1) Use kOT, the energy landscape and kf to determine kPR:
    kPR ,_= calc_rate_PAM_to_Rloop(kf, kOT, energy_landscape)

    # 2) Calculate chi-squared with single mismatch off-targets:
    Epsilon = np.diff(energy_landscape)
    Epsilon[1:] *= -1

    new_parameters = list(Epsilon) + list(epsilon_I)
    # kSP should be irrelevant once equil. between P and S is reached

    # new_parameters.append(np.log10(10 ** 3))
    new_parameters.append(np.log10(0.25))
    new_parameters.append(np.log10(kPR))
    new_parameters.append(np.log10(kf))
    new_parameters = np.array(new_parameters)

    ONtarget_occupancy = 1.0

    V = 0
    for i in range(len(xdata)):
        V += dCas9.calc_Chi_square(new_parameters, xdata[i], ydata[i], yerr[i], ONtarget_occupancy,
                              model_id=model_id)

    return V, new_parameters


def optimize_forward_rate(simset, forward_rates, mode='mean'):
    xdata, ydata, yerr = Bdata.prepare_multiprocessing(use_single_mm_only=True,
                                                       use_on_rate=True,
                                                       use_off_rate=False,
                                                       use_occupancy=False)

    if mode == 'mean':
        landscape_avg, epsI_avg, kOT = average_solution(simset)
    elif mode =='median':
        landscape_avg, epsI_avg, kOT = median_solution(simset)

    # kOT = calc_on_rate()

    _, lower_bnd_kf = calc_rate_PAM_to_Rloop(kf=1.0, kOT=kOT, energy_landscape=landscape_avg)

    V = []
    kf_vals = []
    for kf in forward_rates:
        if kf < lower_bnd_kf:
            continue
        kf_vals.append(kf)
        v, _ = test_forward_rate(kf=kf,
                                 kOT=kOT,
                                 energy_landscape=landscape_avg,
                                 epsilon_I=epsI_avg,
                                 xdata=xdata,
                                 ydata=ydata,
                                 yerr=yerr)
        V.append(v)


    # Now find the optimum:
    kf_opt = kf_vals[np.argmin(V)]
    V_opt, parameters_opt = test_forward_rate(kf=kf_opt,
                                      kOT=kOT,
                                      energy_landscape=landscape_avg,
                                      epsilon_I=epsI_avg,
                                      xdata=xdata,
                                      ydata=ydata,
                                      yerr=yerr)

    return V, kf_vals, kf_opt, parameters_opt, V_opt



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
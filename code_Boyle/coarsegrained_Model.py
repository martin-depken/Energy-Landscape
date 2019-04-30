import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import sys
import CRISPR_dCas9_binding_curve_Boyle as dCas9
sys.path.append('../code_Pclv/')
import CRISPR_Kinetic_model as Pclv
import Weighted_Average as WA

def calc_coarsegrained_Model(parameters, mismatch_positions, model_id='general_energies_no_kPR'):
    t = np.array([500, 1000, 1500])
    epsilon, kf = dCas9.unpack_parameters(parameters, model_id=model_id)
    energies = dCas9.get_energies(epsilon, mismatch_positions)
    kb = dCas9.get_backward_rates(energies, kf)
    kf[-1] = 10 ** 5  ## Adds an absorbing state in the end of the Rloop.
    e_pam = energies[0] + np.log(10)
    energies[0] = e_pam  # convert to 1nM
    ksp = kf[0] / 10  # convert to 1nM
    kf[0] = ksp
    E = np.cumsum(energies)
    Rloop_start = 10
    F = -np.log(np.sum(np.exp(-E[Rloop_start:])))
    Rloop_end = int(np.dot(np.arange(Rloop_start, 21), np.exp(-E[Rloop_start:])) / np.exp(-F))
    Delta = np.log(kf[1:Rloop_end + 2] / kb[1:Rloop_end + 2])
    Delta[-1] = 10 ## Adds an absorbing state after the bound state to calculate the first-passage rate k_R
    k_R = ksp * Pclv.Pclv(Delta)
    k_UB = k_R * np.exp(F)
    k = k_UB + k_R
    P_eq = k_R / k
    P_R = P_eq * (1 - np.exp(-k * t))
    P = np.exp(-e_pam) / (1 + np.exp(-e_pam)) * (1 - P_R) + P_R
    k_on = dCas9.least_squares_line_through_origin(x_points=t, y_points=P)

    return k_on, k_R, k_UB


def predict_coarsegrained_Model(parameters, model_id='general_energies_no_kPR'):
    Ng = 20
    k_on_mat = np.zeros((Ng, Ng))
    k_R_mat = np.zeros((Ng, Ng))
    k_UB_mat = np.zeros((Ng, Ng))

    for i in range(1, Ng + 1):
        for j in range(1, Ng + 1):
            mismatch_positions = [i, j]
            if i == j:
                mismatch_positions = [i]
            k_on, k_R, k_UB = calc_coarsegrained_Model(parameters, mismatch_positions, model_id=model_id)
            k_on_mat[20 - i, 20 - j] = k_on
            k_R_mat[20 - i, 20 - j] = k_R
            k_UB_mat[20 - i, 20 - j] = k_UB

    return k_on_mat, k_R_mat, k_UB_mat

def correlation_coarsegrained_Model(parameters, model_id='general_energies_no_kPR', path='../Data_Boyle/', replica='1', Plot=True):
    wa = WA.calc_Weighted_average(path='../Data_Boyle/', replica='1', save=False)
    prediction = wa[['MM_pos', 'WA_kon']].copy()
    prediction['model_kon'] = wa['MM_pos'].apply(lambda x: calc_coarsegrained_Model(parameters,x,model_id=model_id)[0])
    score = prediction.dropna().apply(lambda x: np.abs(x['WA_kon'] - x['model_kon']) / x['WA_kon'], axis=1).mean()
    corr = prediction[['WA_kon', 'model_kon']].corr().WA_kon[1]
    if Plot:
        plt.figure()
        ymax = np.nanmax([np.nanmax([prediction.WA_kon * 10 ** 4]), np.nanmax([prediction.model_kon * 10 ** 4])])
        plt.plot(prediction.WA_kon * 10 ** 4, prediction.model_kon * 10 ** 4, 'ro')
        plt.plot([0.0, ymax], [0.0, ymax], 'k', lw=2)
        plt.title('Associaton rate 1nM ($10^{-4} s^{-1}$)\n Correlation='+str(float(round(1000*corr))/1000), fontsize=15)
        plt.xlabel('Training data (weighted average)', fontsize=15)
        plt.ylabel('Coarse grained Model', fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)

    return score, corr, prediction

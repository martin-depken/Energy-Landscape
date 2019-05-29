import pandas as pd
import numpy as np
import sys
sys.path.append('../predict_pclv_Stijn/')
import Nucleaseq_data_processing as process
import calculate_cleavage_rate as CLV
sys.path.append('../code_Boyle/')
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style('ticks')
current_colors = sns.color_palette()
sns.set_palette('Accent')

'''
Helper Functions
'''


def Weighted_average(row):
    y = np.array(row['Log_clv'])
    e = np.array(row['error'])
    wa = np.nan
    if len(y) > 0:
        wa = np.average(y, weights=e ** -2, axis=0)
    return wa


def calc_log_kclv_from_Boyle(mismatch_positions, Boyle_param, clv_rate):
    parameters = Boyle_param[1:].copy()
    parameters = np.delete(parameters, -2)
    parameters = np.append(parameters, np.log10(clv_rate))
    model_id = 'Clv_Saturated_general_energies_v2'
    kclv = CLV.calc_clv_rate_fast(parameters, model_id, mismatch_positions)
    return np.log10(kclv)

'''
Main
'''

def predict_clv(Boyle_param, clv_rate=1000.0, Plot=True):
    ## Watch out: for now this only works with
    ##the Boyle parameters which are related to
    ##the model_id='general_energies_no_kPR'

    filename = 'ECas9_cleavage_rate_and_y0_Canonical_OT-r_0-2.csv'
    path = '../data_nucleaseq_Finkelsteinlab/targetE/'
    x, y, e = process.prepare_multiprocessing_nucleaseq_log(filename, path)
    data = pd.DataFrame(data={'MM_pos': x, 'Log_clv': y, 'error': e})
    wa = []
    for i in data.index:
        wa.append(Weighted_average(data.loc[i]))
    data['WA'] = wa
    data['prediction'] = data['MM_pos'].apply(lambda x: calc_log_kclv_from_Boyle(x, Boyle_param, clv_rate))
    score = ((data.prediction - data.WA).apply(np.abs) / data.WA.apply(np.abs)).mean()
    corr = data[['WA', 'prediction']].corr()['WA'].iloc[1]
    chi_square = 0
    for i in data.index:
        y = np.array(data.Log_clv.iloc[i])
        pred = np.array(data.prediction.iloc[i])
        e = np.array(data.error.iloc[i])
        chi_square += np.sum(((y - pred) / e) ** 2)
    data['MM_ID'] = data.MM_pos.apply(lambda x: '|'.join(map(str, x)))
    data.set_index('MM_ID', inplace=True)
    single_dat = data[data.MM_pos.apply(lambda x: len(x) == 1)][['WA', 'prediction']].reset_index(drop=False)
    single_dat['MM_ID'] = pd.to_numeric(single_dat['MM_ID'], downcast='integer')
    single_dat.sort_values('MM_ID', inplace=True)
    double_mat = np.empty((20, 20))
    double_mat[:] = np.nan
    for i in range(20):
        for j in range(i + 1, 20):
            ID = '|'.join(map(str, [i + 1, j + 1]))
            double_mat[i, j] = data.WA.loc[ID]
            double_mat[j, i] = data.prediction.loc[ID]
    if Plot:
        plt.figure()
        ymax = np.nanmax([np.nanmax(data.WA), np.nanmax(data.prediction)])
        ymin = np.nanmin([np.nanmin(data.WA), np.nanmin(data.prediction)])
        plt.plot(data.WA, data.prediction, 'ro')
        plt.plot([ymin, ymax], [ymin, ymax], 'k', lw=2)
        plt.title('log10 Cleavage rate ($s^{-1}$)\n Correlation=' + str(float(round(100 * corr)) / 100), fontsize=15)
        plt.xlabel('data', fontsize=15)
        plt.ylabel('Model prediction', fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        ax = plt.gca()
        ax.set_aspect('equal')
        plt.xlim([ymin, ymax])
        plt.ylim([ymin, ymax])
        plt.figure()
        ax = sns.heatmap(double_mat, cmap='coolwarm', vmax=data.WA.max(), vmin=data.WA.min())
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        labels = ['', '2', '', '4', '', '6', '', '8', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20']
        ax.set_xticklabels(labels);
        ax.set_yticklabels(labels);
        cax = plt.gcf().axes[-1]
        cax.tick_params(labelsize=15)
        plt.title('log10 Cleavage rate ($s^{-1}$)\n Top: Data - Bottom: Model', fontsize=15);
        plt.figure()
        plt.plot(single_dat.MM_ID, single_dat.WA, marker='o', label='Data')
        plt.plot(single_dat.MM_ID, single_dat.prediction, marker='o', label='Model')
        ax = plt.gca()
        ax.set_xlim(1, 21)
        ax.set_xlabel('mismatch position', fontsize=15)
        ax.set_ylabel('log10 Cleavage rate ($s^{-1}$)', fontsize=15)
        ax.set_xticks([i + 0.5 for i in range(20)])
        ax.set_xticklabels(
            [1, '', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0, fontsize=15);
        ax.set_yticklabels(ax.get_yticks(), fontsize=15)
        ax.legend(loc='best', fontsize=12, frameon=True)
        plt.title('Single Mismatches', fontsize=15)

    return chi_square, corr, score, data, single_dat, double_mat


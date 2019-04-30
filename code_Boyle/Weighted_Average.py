import os
import numpy as np
import pandas as pd
import Boyle_data_processing as process
import CRISPR_dCas9_binding_curve_Boyle as dCas9
import matplotlib.pylab as plt


def Weighted_average(row):
    y = row['ydata']
    e = row['yerr']
    wa = np.zeros(3)
    for n in range(3):
        wa[n] = np.nan
        if len(y[n]) > 0:
            wa[n] = np.average(y[n], weights=e[n] ** -2, axis=0)
    return wa


def calc_Weighted_average(path='../Data_Boyle/', replica='1', outputdirectory='../Data_Boyle/Weighted_Avarage/',
                          save=True):
    xdata, ydata, yerr = process.prepare_multiprocessing(replica=replica, path=path)
    data = pd.DataFrame(columns=['xdata', 'ydata', 'yerr'])
    data['xdata'] = xdata
    data['ydata'] = ydata
    data['yerr'] = yerr
    data['WA'] = data.apply(Weighted_average, axis=1)

    WA_data = pd.DataFrame(columns=['MM_pos', 'WA_kon', 'WA_occ', 'WA_koff'])
    WA_data.MM_pos = xdata
    WA_data.WA_kon = data.WA.apply(lambda x: x[1])
    WA_data.WA_occ = data.WA.apply(lambda x: x[0])
    WA_data.WA_koff = data.WA.apply(lambda x: x[2])

    occ = np.zeros((20, 20))
    kon = np.zeros((20, 20))
    koff = np.zeros((20, 20))
    for ind in data.index:
        x = data['xdata'].iloc[ind]
        WA = data['WA'].iloc[ind]
        if len(x) > 0:
            n1 = 20 - x[0]
            n2 = 20 - x[0]
            if len(x) > 1:
                n2 = 20-x[1]
            occ[n1, n2] = WA[0]
            occ[n2, n1] = WA[0]
            kon[n1, n2] = WA[1]
            kon[n2, n1] = WA[1]
            koff[n1, n2] = WA[2]
            koff[n2, n1] = WA[2]
        else:
            occ_OT = WA[0]
            kon_OT = WA[1]
            koff_OT = WA[2]
    if save:
        if not os.path.exists(outputdirectory):
            os.makedirs(outputdirectory)
        np.savetxt(outputdirectory + 'DataOccupancy.txt', occ, delimiter=',')
        np.savetxt(outputdirectory + 'DataOnRate.txt', kon, delimiter=',')
        np.savetxt(outputdirectory + 'DataOffRate.txt', koff, delimiter=',')
        np.savetxt(outputdirectory + 'DataOccupancy_OT.txt', np.array([occ_OT]), delimiter=',')
        np.savetxt(outputdirectory + 'DataOnRate_OT.txt', np.array([kon_OT]), delimiter=',')
        np.savetxt(outputdirectory + 'DataOffRate_OT.txt', np.array([koff_OT]), delimiter=',')

    return WA_data


def predict_train(parameters, model_id='general_energies_no_kPR', path='../Data_Boyle/', replica='1', Plot=True):
    wa = calc_Weighted_average(path='../Data_Boyle/', replica='1', save=False)
    prediction = wa[['MM_pos', 'WA_kon']].copy()
    prediction['model_kon'] = wa['MM_pos'].apply(lambda x: dCas9.calc_Boyle(False, False, True, parameters, x, model_id=model_id)[1])
    score = prediction.dropna().apply(lambda x: np.abs(x['WA_kon'] - x['model_kon']) / x['WA_kon'], axis=1).mean()
    corr = prediction[['WA_kon', 'model_kon']].corr().WA_kon[1]
    if Plot:
        plt.figure()
        ymax = np.nanmax([np.nanmax([prediction.WA_kon * 10 ** 4]), np.nanmax([prediction.model_kon * 10 ** 4])])
        plt.plot(prediction.WA_kon * 10 ** 4, prediction.model_kon * 10 ** 4, 'ro')
        plt.plot([0.0, ymax], [0.0, ymax], 'k', lw=2)
        plt.title('Associaton rate 1nM ($10^{-4} s^{-1}$)\n Correlation='+str(float(round(1000*corr))/1000), fontsize=15)
        plt.xlabel('Training data (weighted average)', fontsize=15)
        plt.ylabel('Model prediction', fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)

    return score, corr, prediction

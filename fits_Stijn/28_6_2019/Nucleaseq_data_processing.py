import pandas as pd
import numpy as np
import sys
sys.path.append('../code_Boyle/')
sys.path.append('/home/svandersmagt/Energy_Landscape_dCas9/code_general/')
import Boyle_data_processing as processing

def prepare_multiprocessing_nucleaseq(filename, path, fit_to_median=False):
    data = pd.read_csv(path + filename,
                       usecols=['Mutation Positions', 'Mutation ID',
                                'cleavage_rate', 'cleavage_rate_5th_pctl', 'cleavage_rate_95th_pctl', ],
                       na_filter=False)
    data['error'] = data.apply(lambda x: np.max([x['cleavage_rate'] - x['cleavage_rate_5th_pctl'],
                                                 x['cleavage_rate_95th_pctl'] - x['cleavage_rate']]), axis=1)
    data = data[['cleavage_rate', 'error', 'Mutation Positions']]
    if fit_to_median:
        data['sq_error'] = data['error'] ** 2
        grouped_data = pd.DataFrame(columns=['Mutation Positions', 'cleavage_rate'])
        grouped_data = data.groupby('Mutation Positions')[
            ['Mutation Positions', 'cleavage_rate']].median().reset_index()
        grouped_data['sum_sq_error'] = data.groupby('Mutation Positions')['sq_error'].sum().reset_index()['sq_error']
        grouped_data['seq_num'] = data.groupby('Mutation Positions')['sq_error'].agg(lambda x: len(x)).reset_index()[
            'sq_error']
        grouped_data['sq_mean_error'] = grouped_data['sum_sq_error'] / grouped_data['seq_num'] ** 2
        grouped_data['error'] = grouped_data['sq_mean_error'].apply(lambda x: np.sqrt(x))
        grouped_data = grouped_data[['cleavage_rate', 'error', 'Mutation Positions']]
        grouped_data['error'] = grouped_data['error'].apply(lambda x: [x])
        grouped_data['cleavage_rate'] = grouped_data['cleavage_rate'].apply(lambda x: [x])
    else:
        grouped_data = data.groupby('Mutation Positions').agg(lambda x: list(x)).reset_index()
    grouped_data['Mutation Positions list'] = grouped_data['Mutation Positions'].apply(
        lambda x: map(int, x.split('|')) if not x == '' else [])
    xdata = grouped_data['Mutation Positions list'].tolist()
    ydata = grouped_data['cleavage_rate'].tolist()
    yerr = grouped_data['error'].tolist()

    return xdata, ydata, yerr

def weighting(yerr): #used for calculating weighted average
    yerr_sqr = np.zeros(len(yerr))
    for i in range(len(yerr)):
        yerr_sqr[i] = yerr[i]**2
    Z = np.sum(np.reciprocal(yerr_sqr))
    weights = np.zeros(len(yerr))
    for i in range(len(yerr)):
        weights[i] = 1/yerr_sqr[i]/Z
    error_of_wa = np.sqrt(1/Z)
    return weights, error_of_wa

def prepare_multiprocessing_combined(rep_on,filename_clv,path_on,path_clv,fit_to_wa=False,only_single=False,log_boyle=False,fill_in_boyle=False):
    
    reload(processing)
    xdata_clv, ydata_clv, yerr_clv = prepare_multiprocessing_nucleaseq_log(filename_clv,path_clv)
    xdata_on, ydata_on, yerr_on = processing.prepare_multiprocessing(rep_on,path_on,True,False,False,False,False)
    
    ydata = list()
    yerr = list()
    for i in range(len(xdata_clv)):
        for j in range(len(xdata_on)):
            if (len(xdata_clv[i])==len(xdata_on[j])):
                if(len(xdata_clv[i])==0 and len(xdata_on[j])==0):
                    ydata.append([ydata_clv[i],ydata_on[j][1]])
                    yerr.append([yerr_clv[i],yerr_on[j][1]])
                if(len(xdata_clv[i])==1 and xdata_clv[i][0]==xdata_on[j][0]):
                    ydata.append([ydata_clv[i],ydata_on[j][1]])
                    yerr.append([yerr_clv[i],yerr_on[j][1]])
                if((len(xdata_clv[i])==2 and xdata_clv[i][0]==xdata_on[j][0] and xdata_clv[i][1]==xdata_on[j][1])
                  or (len(xdata_clv[i])==2 and xdata_clv[i][0]==xdata_on[j][1] and xdata_clv[i][1]==xdata_on[j][0])):
                    ydata.append([ydata_clv[i],ydata_on[j][1]])
                    yerr.append([yerr_clv[i],yerr_on[j][1]]) 
    if fit_to_wa:
        for i in range(len(xdata_clv)):
            weightsclv, errorclv = weighting(yerr[i][0])
            ydata[i][0] = [np.average(ydata[i][0],weights=weightsclv)]
            if len(ydata[i][1])==0:
                ydata[i][1] = []
            else:
                weightson, erroron = weighting(yerr[i][1])
                ydata[i][1] = [np.average(ydata[i][1],weights=weightson)]
            
            yerr[i][0] = [errorclv]
            if len(ydata[i][1])==0:
                yerr[i][1] = []
            else:
                yerr[i][1] = [erroron]
                
    if only_single:
        xdata_single = []
        ydata_single = []
        yerr_single = []
        for i in range(len(xdata_clv)):
            if len(xdata_clv[i])<2:
                xdata_single.append(xdata_clv[i])
                ydata_single.append(ydata[i])
                yerr_single.append(yerr[i])
        return xdata_single, ydata_single, yerr_single
    
    if log_boyle:
        for i in range(len(xdata_clv)):
            for j in range(len(ydata[i][1])):
                datapoint = ydata[i][1][j]
                errorpoint = yerr[i][1][j]
                ydata[i][1][j] = np.log10(datapoint)
                yerr[i][1][j] = np.log10(1 + 2*errorpoint/(datapoint-2*errorpoint))
                
    if fill_in_boyle and  fit_to_wa:
        for a in range(1,20):
            temp = []
            temperror = []
            for i in range(len(xdata_clv)):
                if len(xdata_clv[i])==2 and xdata_clv[i][0]==a:
                    if not len(ydata[i][1])==0:
                        temp.append(ydata[i][1][0])
                        temperror.append(yerr[i][1][0])
            
            avg = np.average(temp,weights=np.reciprocal(temperror))
            err = np.max(temperror)
            
            for i in range(len(xdata_clv)):
                if len(xdata_clv[i])==2 and xdata_clv[i][0]==a:
                    if len(ydata[i][1])==0:
                        ydata[i][1] = [avg]
                        yerr[i][1] = [err]
            
    
    return xdata_clv, ydata, yerr
    
def prepare_multiprocessing_nucleaseq_log(filename, path, fit_to_wa=False):
    data = pd.read_csv(path + filename,
                       usecols=['Mutation Positions', 'Mutation ID',
                                'cleavage_rate', 'cleavage_rate_5th_pctl', 'cleavage_rate_95th_pctl'],
                       na_filter=False)
    data['cleavage_rate_log'] = data['cleavage_rate'].apply(np.log10)
    data['cleavage_rate_5th_pctl_log'] = data['cleavage_rate_5th_pctl'].apply(np.log10)
    data['cleavage_rate_95th_pctl_log'] = data['cleavage_rate_95th_pctl'].apply(np.log10)
    data['error_log'] = data.apply(lambda x: np.max([x['cleavage_rate_log'] - x['cleavage_rate_5th_pctl_log'],
                                                     x['cleavage_rate_95th_pctl_log'] - x['cleavage_rate_log']]),
                                   axis=1)
    data = data[['cleavage_rate_log', 'error_log', 'Mutation Positions']]
    grouped_data = data.groupby('Mutation Positions').agg(lambda x: list(x)).reset_index()
    grouped_data['Mutation Positions list'] = grouped_data['Mutation Positions'].apply(
        lambda x: map(int, x.split('|')) if not x == '' else [])
    xdata = grouped_data['Mutation Positions list'].tolist()
    ydata = grouped_data['cleavage_rate_log'].tolist()
    yerr = grouped_data['error_log'].tolist()
    
    if fit_to_wa:
        for i in range(len(xdata)):
            weightsclv, errorclv = weighting(yerr[i])
            ydata[i] = [np.average(ydata[i],weights=weightsclv)]
            yerr[i] = [errorclv]
        
    return xdata, ydata, yerr

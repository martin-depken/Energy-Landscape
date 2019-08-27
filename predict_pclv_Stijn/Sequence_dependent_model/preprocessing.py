import pandas as pd
import numpy as np
import sys
sys.path.append('../code_Boyle/')
sys.path.append('/home/svandersmagt/Energy_Landscape_dCas9/code_general/')



def prepare_multiprocessing_seq_dep(filename, path):
    data = pd.read_csv(path + filename,
                       usecols=['Sequence','Target','Mutation Positions', 'Mutation ID',
                                'cleavage_rate', 'cleavage_rate_5th_pctl', 'cleavage_rate_95th_pctl'],
                       na_filter=False)
    data['cleavage_rate_log'] = data['cleavage_rate'].apply(np.log10)
    data['cleavage_rate_5th_pctl_log'] = data['cleavage_rate_5th_pctl'].apply(np.log10)
    data['cleavage_rate_95th_pctl_log'] = data['cleavage_rate_95th_pctl'].apply(np.log10)
    data['error_log'] = data.apply(lambda x: np.max([x['cleavage_rate_log'] - x['cleavage_rate_5th_pctl_log'],
                                                     x['cleavage_rate_95th_pctl_log'] - x['cleavage_rate_log']]),
                                   axis=1)
    data = data[['cleavage_rate_log', 'error_log', 'Mutation Positions','Sequence','Target']]
    grouped_data = data.groupby('Sequence').agg(lambda x: list(x)).reset_index()
    grouped_data['Target'] = grouped_data['Target'].apply(lambda x: x[0])
    grouped_data['Mutation Positions list'] = grouped_data['Mutation Positions'].apply(
        lambda x: (map(int, (str(x).split('|'))  if not x == '' else []) if not type(x) == type(list()) else (map(int, str(x[0]).split('|')) if not x[0] == '' else [])))
    grouped_data['Sequence_Target'] = grouped_data.apply(
        lambda s: [s['Sequence'],s['Target'],s['Mutation Positions list']],axis=1)
    xdata = grouped_data['Sequence_Target'].tolist()
    ydata = grouped_data['cleavage_rate_log'].tolist()
    yerr = grouped_data['error_log'].tolist()
    return xdata, ydata, yerr

def prepare_multiprocessing_seq_dep_aba(filename, path):
    data = pd.read_csv(path + filename,
                       usecols=['sequence','Target','Mutation Positions',
                                'ABA', 'error'],
                       na_filter=False)
    grouped_data = data.groupby('sequence').agg(lambda x: list(x)).reset_index()
    grouped_data['Target'] = grouped_data['Target'].apply(lambda x: x[0])
    grouped_data['Mutation Positions list'] = grouped_data['Mutation Positions'].apply(
        lambda x: (map(int, (str(x).split('|'))  if not x == '' else []) if not type(x) == type(list()) else (map(int, str(x[0]).split('|')) if not x[0] == '' else [])))
    grouped_data['Sequence_Target'] = grouped_data.apply(
        lambda s: [s['sequence'],s['Target'],s['Mutation Positions list']],axis=1)
    xdata = grouped_data['Sequence_Target'].tolist()
    ydata = grouped_data['ABA'].tolist()
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

def prepare_multiprocessing_seq_dep_wa(filename, path):
    xdata,ydata_full,yerr_full = prepare_multiprocessing_seq_dep(filename, path)
    ydata = []
    yerr = []
    for i in range(len(xdata)):
        weights,error_of_wa = weighting(yerr_full[i])
        ydata.append([np.average(ydata_full[i],weights=weights)])
        yerr.append([error_of_wa])
    return xdata, ydata, yerr

def prepare_multiprocessing_seq_dep_wa_aba(filename, path):
    xdata,ydata_full,yerr_full = prepare_multiprocessing_seq_dep_aba(filename, path)
    ydata = []
    yerr = []
    for i in range(len(xdata)):
        weights,error_of_wa = weighting(yerr_full[i])
        ydata.append([np.average(ydata_full[i],weights=weights)])
        yerr.append([error_of_wa])
    return xdata, ydata, yerr

def prepare_multiprocessing_seq_indep(filename, path):
    data = pd.read_csv(path + filename,
                       usecols=['Mutation Positions','cleavage_rate', 'cleavage_rate_5th_pctl', 'cleavage_rate_95th_pctl'],
                       na_filter=False)
    data['cleavage_rate_log'] = data['cleavage_rate'].apply(np.log10)
    data['cleavage_rate_5th_pctl_log'] = data['cleavage_rate_5th_pctl'].apply(np.log10)
    data['cleavage_rate_95th_pctl_log'] = data['cleavage_rate_95th_pctl'].apply(np.log10)
    data['error_log'] = data.apply(lambda x: np.max([x['cleavage_rate_log'] - x['cleavage_rate_5th_pctl_log'],
                                                     x['cleavage_rate_95th_pctl_log'] - x['cleavage_rate_log']]),
                                   axis=1)
    data = data[['cleavage_rate_log', 'error_log', 'Mutation Positions']]
    grouped_data = data.groupby('Mutation Positions').agg(lambda x: list(x)).reset_index()
    grouped_data['Mutation Positions list'] = grouped_data['Mutation Positions'].apply(lambda x: (map(int, (x.split('|'))  if not x == '' else []) if not type(x) == type(list()) else (map(int, x[0].split('|')) if not x[0] == '' else [])))
    xdata = grouped_data['Mutation Positions list'].tolist()
    ydata = grouped_data['cleavage_rate_log'].tolist()
    yerr = grouped_data['error_log'].tolist()
    return xdata, ydata,yerr
    
                                                        
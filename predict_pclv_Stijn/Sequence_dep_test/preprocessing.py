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
    grouped_data = data.groupby(['Sequence','Target']).agg(lambda x: list(x)).reset_index()
    grouped_data['Mutation Positions list'] = grouped_data['Mutation Positions'].apply(
        lambda x: (map(int, (x.split('|'))  if not x == '' else []) if not type(x) == type(list()) else (map(int, x[0].split('|')) if not x[0] == '' else [])))
    grouped_data['Sequence_Target'] = grouped_data.apply(
        lambda s: [s['Sequence'],s['Target'],s['Mutation Positions list']],axis=1)
    xdata = grouped_data['Sequence_Target'].tolist()
    ydata = grouped_data['cleavage_rate_log'].tolist()
    yerr = grouped_data['error_log'].tolist()
    return xdata, ydata, yerr
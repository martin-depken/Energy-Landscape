import pandas as pd
import numpy as np

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

def prepare_multiprocessing_nucleaseq_log(filename, path):
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
    return xdata, ydata, yerr

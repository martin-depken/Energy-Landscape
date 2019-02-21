import pandas as pd
import numpy as np

def prepare_multiprocessing_nucleaseq(filename, path):
    data = pd.read_csv(path + filename,
                       usecols=['Mutation Positions', 'Mutation ID',
                                'cleavage_rate', 'cleavage_rate_5th_pctl', 'cleavage_rate_95th_pctl', ],
                       na_filter=False)
    data['error'] = data.apply(lambda x: np.max([x['cleavage_rate'] - x['cleavage_rate_5th_pctl'],
                                                 x['cleavage_rate_95th_pctl'] - x['cleavage_rate']]), axis=1)
    data = data[['cleavage_rate', 'error', 'Mutation Positions']]
    grouped_data = data.groupby('Mutation Positions').agg(lambda x: list(x)).reset_index()
    grouped_data['Mutation Positions list'] = grouped_data['Mutation Positions'].apply(
        lambda x: map(int, x.split('|')) if not x == '' else [])
    grouped_data

    xdata = grouped_data['Mutation Positions list'].tolist()
    ydata = grouped_data['cleavage_rate'].tolist()
    yerr = grouped_data['error'].tolist()

    return xdata, ydata, yerr

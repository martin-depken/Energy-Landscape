import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('ticks')
current_colors = sns.color_palette()

def plot_single_mut_data(data, data_colname, error_colname, Mut_type, Canonical=True, data_name='Data', Plot=True, logplot=False,
                         SaveFigures=False, figure_name='./Figure.pdf'):
    '''
    Plots the measured quantity for sequences containing a single mutation versus the position of the mutation
    (median is taken over the sequence which corresponds to the mutation).

    :param data: The processed data (from Process_SeqLibrary_Finkelsteinlab) loaded in a pandas data frame.
    :param data_colname: The name of the column in the data to be plotted.
    :param error_colname: The name of the column in the data which contains the errors.
    :param Mut_type: Type of the mutation ['r','i', 'd'].
    :param Canonical: Boolian; If True, keeps the mutant with Canonical PAM, if False, keeps only the mutants with
                      noncanonical PAM.
    :param data_name: The name you want to appear on the y-axis and on the Fig. title.
    :param Plot: Boolian; whether or not to produce a plot.
    :param logplot: Boolian; whether or not set the y-axis on log scale.
    :param SaveFigures: Boolian; whether or not to save the Fig.
    :param figure_name: Path+filename with which the file is being saved.
    :return:
            single_mut_data_median: A data frame contains the measurements for all sequence with a single mutation
            (median is taken over the sequence).
    '''

    Mut_names = {'r': 'Mismatch', 'i': 'Insertion', 'd': 'Deletion'}
    Mut_name = Mut_names[Mut_type]
    if error_colname is not None:
        cols = ['Mutation Positions','Mutation Count',data_colname, error_colname]
    else:
        cols = ['Mutation Positions', 'Mutation Count', data_colname]
    single_mut_data = data[(data['Canonical'] == Canonical) & (data['Mutation Type'] == Mut_type) &
                           (data['Mutation Count'] == 1)][cols].copy()
    single_mut_data_WA = Calc_Weighted_average(single_mut_data, data_colname, error_colname)
    single_mut_data_WA['Position'] = single_mut_data_WA['Mutation Positions'].apply(lambda x: get_mut_pos(x)[0])
    single_mut_data_WA.sort_values(by='Position', inplace=True)
    single_mut_data_WA.reset_index(drop=True, inplace=True)
    single_mut_data_WA.drop('Mutation Positions', axis=1, inplace=True)
    wa_colname = data_colname + ' (weighted average)'
    if Plot:
        if logplot:
            plt.semilogy(single_mut_data_WA['Position'], single_mut_data_WA[wa_colname]
                         , marker='o', linewidth=1, color='blue')
        else:
            plt.plot(single_mut_data_WA['Position'], single_mut_data_WA[wa_colname]
                         , marker='o', linewidth=1, color='blue')
        plt.xlabel(Mut_name + ' postion', fontsize=15)
        plt.ylabel(data_name, fontsize=15)
        plt.title(data_name, fontsize=15)
        plt.legend(loc='best')
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        if not Canonical:
            plt.title(data_name + ' - Noncanonical', fontsize=15)
        if SaveFigures:
            plt.savefig(figure_name,format='pdf',bbox_inches='tight')

    return single_mut_data_WA


def plot_double_mut_data(data, data_colname, error_colname, Mut_type, Canonical=True, Ng=20, data_name='Data', Plot=True,
                        logplot=False, SaveFigures=False, figure_name='./Figure.pdf'):

    '''
    Plots the heat-map showing the musured quantity for all sequences containing two mutations (of the same type).
    (median is taken over the sequence which corresponds to the mutation).

    :param data:  The processed data (from Process_SeqLibrary_Finkelsteinlab) loaded in a pandas data frame.
    :param data_colname: The name of the column in the data to be plotted.
    :param error_colname: The name of the column in the data which contains the errors.
    :param Mut_type: Type of the mutation ['r','i', 'd'].
    :param Canonical: Boolian; If True, keeps the mutant with Canonical PAM, if False, keeps only the mutants with
                      noncanonical PAM.
    :param Ng: Length of the guid-target hybrid (excluding PAM).
    :param data_name: The name you want to appear on the y-axis and on the Fig. title.
    :param Plot: Boolian; whether or not to produce a plot.
    :param logplot: Boolian; if True, log10 of the measured quantity is plotted.
    :param SaveFigures: Boolian; whether or not to save the Fig.
    :param figure_name: Path+filename with which the file is being saved.
    :return:
            double_mut_data_median: A data frame containing the measurements for all sequence with double mutations
            (median is taken over the sequence).

            double_mut_df: The same as double_mut_data_median, but in the matrix form: each column and row of the
            data frame corresponds to a mismatch position. This is ready to be passed to the sns.heatmap.
            NOTE: if logplot = True, the log10 of the measured quantity is used to build this data frame.
    '''

    Mut_names = {'r': 'Mismatch', 'i': 'Insertion', 'd': 'Deletion'}
    Mut_name = Mut_names[Mut_type]
    if error_colname is not None:
        cols = ['Mutation Positions','Mutation Count',data_colname, error_colname]
    else:
        cols = ['Mutation Positions', 'Mutation Count', data_colname]
    double_mut_data = data[(data['Canonical'] == Canonical) & (data['Mutation Type'] == Mut_type) &
                           (data['Mutation Count'] == 2)][cols].copy()
    double_mut_data_WA = Calc_Weighted_average(double_mut_data, data_colname, error_colname)
    wa_colname = data_colname + ' (weighted average)'
    double_mut_map = np.empty((Ng, Ng))
    double_mut_map[:] = np.nan
    for n in range(len(double_mut_data_WA)):
        pos = double_mut_data_WA['Mutation Positions'].iloc[n]
        Positions = map(lambda x: int(x) - 1, pos.split('|'))
        datapoint = double_mut_data_WA[wa_colname].iloc[n]
        if logplot:
            datapoint = np.log10(datapoint)
        double_mut_map[Positions[1], Positions[0]] = datapoint
    double_mut_df = pd.DataFrame(double_mut_map, index=range(1, Ng+1), columns=range(1, Ng+1))

    if Plot:
        plt.figure()
        sns.heatmap(double_mut_df, cmap='coolwarm', cbar=True)
        plt.xlabel(Mut_name + ' 1', fontsize=15)
        plt.ylabel(Mut_name + ' 2', fontsize=15)
        title = data_name
        if logplot:
            title = 'log10 (' + title + ')'
        plt.title(title, fontsize=15)
        if not Canonical:
            plt.title(title + ' - Noncanonical', fontsize=15)
        if SaveFigures:
            plt.savefig(figure_name,format='pdf',bbox_inches='tight')
            
    return double_mut_data_WA, double_mut_df

def plot_block_mm_data(data, data_colname, error_colname, data_name='Data', Canonical=True, Ng=20, Plot=True, logplot=False,
                        SaveFigures=False, figure_names=['./Figure1.pdf','./Figure2.pdf']):
    '''
    Plots the measured quantity for all mutants containing a mismatch block
    (median is taken over the sequence of the block).
    Also plots a 1D graph where median is taken over the block ending position.

    :param data: The processed data (from Process_SeqLibrary_Finkelsteinlab) loaded in a pandas data frame.
    :param data_colname: The name of the column in the data to be plotted.
    :param error_colname: The name of the column in the data which contains the errors.
    :param data_name: The name you want to appear on the y-axis and on the Fig. title.
    :param Canonical: Boolian; If True, keeps the mutant with Canonical PAM, if False, keeps only the mutants with
                      noncanonical PAM.
    :param Ng: Length of the guid-target hybrid (excluding PAM).
    :param Plot: Boolian; whether or not to produce a plot.
    :param logplot: Boolian; if True, log10 of the measured quantity is plotted, and in the 1D plot the y-axis is set
           to log scale.
    :param SaveFigures: Boolian; whether or not to save the Figs.
    :param figure_names: Path+filenames with which the two Figure files is being saved.
    :return:
            block_mm_data_median: A data frame containing the measurements for all sequences containing mismatch blocks
            (median is taken over the sequence of the blocks).

            block_mm_df: The same as block_mm_data_median, but in the matrix form: each row of the
            data frame corresponds to a block starting position, and each column to the block ending position.
            This is ready to be passed to the sns.heatmap.
            NOTE: if logplot = True, the log10 of the measured quantity is used to build this data frame.

            block_start_mm_data: Build from block_mm_data_median, where median is taken over the block ending position.
    '''
    if error_colname is not None:
        cols = ['Mutation Positions','Mutation Count',data_colname, error_colname]
    else:
        cols = ['Mutation Positions', 'Mutation Count', data_colname]
    multi_mm_data = data[(data['Canonical'] == Canonical) & (data['Mutation Type'] == 'r') &
                           (data['Mutation Count'] > 2)][cols].copy()
    multi_mm_data['All_Positions'] = multi_mm_data['Mutation Positions'].apply(lambda x: get_mut_pos(x))
    multi_mm_data['Position'] = multi_mm_data['All_Positions'].apply(make_block_ID)
    multi_mm_data = multi_mm_data[(multi_mm_data['Position'] != '')]
    if error_colname is not None:
        cols = ['Position', 'Mutation Count', data_colname, error_colname]
    else:
        cols = ['Position', 'Mutation Count', data_colname]
    multi_mm_data = multi_mm_data[cols]
    multi_mm_data.rename(columns={'Position':'Mutation Positions'},inplace=True)
    block_mm_data_WA = Calc_Weighted_average(multi_mm_data, data_colname, error_colname)
    block_mm_data_WA.rename(columns={'Mutation Positions':'Position'},inplace=True)
    wa_colname = data_colname + ' (weighted average)'
    block_mm_map = np.empty((Ng, Ng))
    block_mm_map[:] = np.nan
    for n in range(len(block_mm_data_WA)):
        pos = block_mm_data_WA['Position'].iloc[n]
        Positions = map(lambda x: int(x) - 1, pos.split('|'))
        datapoint = block_mm_data_WA[wa_colname].iloc[n]
        if logplot:
            datapoint = np.log10(datapoint)
        block_mm_map[Positions[0], Positions[1]] = datapoint
    block_mm_df = pd.DataFrame(block_mm_map, index=range(1, Ng + 1), columns=range(1, Ng + 1))

    if Plot:
        plt.figure()
        sns.heatmap(block_mm_df, cmap='coolwarm', cbar=True)
        plt.xlabel('Block end', fontsize=15)
        plt.ylabel('Block start', fontsize=15)
        title = data_name
        if logplot:
            title = 'log10 ('+title+')'
        plt.title(title, fontsize=15)
        if not Canonical:
            plt.title(title+ ' - Noncanonical', fontsize=15)
        if SaveFigures:
            plt.savefig(figure_names[0],format='pdf',bbox_inches='tight')


    block_start_mm_data = block_mm_data_WA.set_index('Position').groupby(
        lambda x: int(x.split('|')[0])).mean().reset_index().rename(columns={'index': 'Block start'})

    if Plot:
        plt.figure()
        if logplot:
            plt.semilogy(block_start_mm_data['Block start'], block_start_mm_data[wa_colname]
                         , marker='o', linewidth=1, color='blue')
        else:
            plt.plot(block_start_mm_data['Block start'], block_start_mm_data[wa_colname]
                         , marker='o', linewidth=1, color='blue')
        plt.xlabel('Block start', fontsize=15)
        plt.ylabel(data_name, fontsize=15)
        plt.title(data_name, fontsize=15)
        if not Canonical:
            plt.title(data_name + ' - Noncanonical', fontsize=15)
        plt.xticks(fontsize=15);
        plt.yticks(fontsize=15);
        if SaveFigures:
            plt.savefig(figure_names[1],format='pdf',bbox_inches='tight')


    return block_mm_data_WA, block_mm_df, block_start_mm_data


'''
Helper functions 
'''
def Weighted_average(row, value_colname, error_colname):
    y = np.array(row[value_colname])
    if error_colname is not None:
        e = np.array(row[error_colname])
    else:
        e = np.ones(len(y))
    wa = np.nan
    if (len(y) > 0):
        wa = np.average(y, weights=e ** -2, axis=0)
    return wa


def Calc_Weighted_average(data, value_colname, error_colname):
    data = data[data['Mutation Count'].apply(lambda x: not np.isnan(x))]
    data.fillna('OT', inplace=True)
    if error_colname is not None:
        data = data[data[error_colname] > 0]
        data = data[['Mutation Positions', value_colname, error_colname]]
        grouped_vals = data[['Mutation Positions', value_colname]].groupby('Mutation Positions').agg(lambda x: list(x))
        grouped_errs = data[['Mutation Positions', error_colname]].groupby('Mutation Positions').agg(lambda x: list(x))
        grouped = grouped_vals.merge(grouped_errs, left_index=True, right_index=True)
    else:
        data = data[['Mutation Positions', value_colname]]
        grouped = data[['Mutation Positions', value_colname]].groupby('Mutation Positions').agg(lambda x: list(x))
    grouped.reset_index(inplace=True)

    wa = []
    for i in grouped.index:
        wa.append(Weighted_average(grouped.loc[i], value_colname, error_colname))
    wa_colname = value_colname + ' (weighted average)'
    grouped[wa_colname] = wa
    if error_colname is not None:
        cols = [value_colname, error_colname]
    else:
        cols = value_colname
    grouped.drop(cols, inplace=True, axis=1)

    return grouped

def make_block_ID(MM_pos_list):
    bi = MM_pos_list[0]
    bf = MM_pos_list[-1]
    block = np.arange(bi, bf + 1)
    if (len(block) != len(MM_pos_list)) or np.any(block != np.array(MM_pos_list)):
        return ''
    return '|'.join(map(str, [bi, bf]))


def get_mut_pos(mut_pos_str):
    if mut_pos_str =='':
        return np.NAN
    return map(int, mut_pos_str.split('|'))


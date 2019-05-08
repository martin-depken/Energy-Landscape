import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('ticks')
current_colors = sns.color_palette()
import importlib as imp
import sys
sys.path.append('../code_general_Finkelsteinlab/')
import plotting_Finkelsteinlab_Diewertje as pltF

import Calculate_ABA_Finkelsteinlab_Diewertje as CalcABA
imp.reload(CalcABA);

def plot_single_mut_ABA(data, Mut_type, data_name='Finkelstein Data', Canonical=True, Plot=True):
    Mut_names = {'r': 'Mismatch', 'i': 'Insertion', 'd': 'Deletion'}
    Mut_name = Mut_names[Mut_type]

    data['Mutation Type'].fillna('', inplace=True)
    single_mut_data = data[(data['Canonical'] == Canonical) & (data['Mutation Type'] == Mut_type)][
        ['Mutation ID', 'Delta ABA (kBT)', 'Uncertainty']]
    single_mut_data['Position'] = single_mut_data['Mutation ID'].apply(lambda x: int(x.split(':')[1]))
    single_mut_data_mean = single_mut_data[['Position', 'Delta ABA (kBT)', 'Uncertainty']].groupby(
        'Position').mean().reset_index()

    if Plot:
        single_mut_data_mean.plot(x='Position', y='Delta ABA (kBT)', yerr='Uncertainty', marker='o', linewidth=1,
                                  color='blue')
        plt.xlabel(Mut_name + ' postion', fontsize=15)
        plt.ylabel('$\Delta$ABA', fontsize=15)
        plt.title(data_name, fontsize=15)
        if not Canonical:
            plt.title(data_name + ' - Noncanonical', fontsize=15)

    return single_mut_data_mean


def plot_double_mut_ABA(data, Mut_type, data_name='Finkelstein Data', Canonical=True, Plot=True):
    Mut_names = {'r': 'Mismatch', 'i': 'Insertion', 'd': 'Deletion'}
    Mut_name = Mut_names[Mut_type]
    Mut_type_str = Mut_type + '|' + Mut_type

    data['Mutation Type'].fillna('', inplace=True)
    double_mut_data = data[(data['Canonical'] == Canonical) & (data['Mutation Type'] == Mut_type_str)][
        ['Mutation ID', 'Delta ABA (kBT)', 'Uncertainty']]
    double_mut_data['Position'] = double_mut_data['Mutation ID'].apply(
        lambda x: '|'.join(map(lambda y: y.split(':')[1], x.split('|'))))
    double_mut_data_mean = double_mut_data[['Position', 'Delta ABA (kBT)', 'Uncertainty']].groupby('Position').mean().reset_index()

    Ng = 20;
    double_mut_map = np.empty((Ng, Ng))
    double_mut_map[:] = np.nan
    for n in range(len(double_mut_data_mean)):
        pos = double_mut_data_mean['Position'].iloc[n]
        Positions = map(lambda x: int(x) - 1, pos.split('|'))
        double_mut_map[Positions[1], Positions[0]] = double_mut_data_mean['Delta ABA (kBT)'].iloc[n]

    if Plot:
        plt.figure()
        sns.heatmap(double_mut_map, cmap='coolwarm', cbar=True, vmin=0,
                    vmax=double_mut_data_mean['Delta ABA (kBT)'].max())
        plt.xlabel(Mut_name + ' 1', fontsize=15)
        plt.ylabel(Mut_name + ' 2', fontsize=15)
        plt.title('$\Delta$ABA - ' + data_name, fontsize=15)
        if not Canonical:
            plt.title('$\Delta$ABA - ' + data_name + ' - Noncanonical', fontsize=15)
        ax = plt.gca()
        ax.set_xticklabels(map(lambda x: str(int(x)), ax.get_xticks() + 0.5));
        ax.set_yticklabels(list(map(lambda x: str(int(x)), ax.get_yticks() + 0.5)));

    return double_mut_data_mean, double_mut_map


def make_block_ID(MM_pos_list):
    bi = MM_pos_list[0]
    bf = MM_pos_list[-1]
    block = np.arange(bi, bf + 1)
    if (len(block) != len(MM_pos_list)) or np.any(block != np.array(MM_pos_list)):
        return ''
    return '|'.join(list(map(str, [bi, bf])))


def plot_block_mm_ABA(data, data_name='Finkelstein Data', Canonical=True, Plot=True,
                      SaveFigures=False, figure_names=['figure1','figure 2']):


    data['Mutation Type'].fillna('', inplace=True)
    select_multi_mm = data['Mutation Type'].apply(lambda x: np.unique(np.array(x.split('|')))[0] == 'r')
    multi_mm_data = data[(data['Canonical'] == Canonical) & (select_multi_mm) & (data['Mutation Count'] > 2)][['Mutation ID', 'Delta ABA (kBT)', 'Uncertainty']]
    multi_mm_data['MM Positions'] = multi_mm_data['Mutation ID'].apply(
        lambda x: list(map(lambda y: int(y.split(':')[1]), x.split('|'))))
    multi_mm_data['Position'] = multi_mm_data['MM Positions'].apply(make_block_ID)
    multi_mm_data = multi_mm_data[(multi_mm_data['Position'] != '')]
    block_mm_data_mean = multi_mm_data[['Position', 'Delta ABA (kBT)', 'Uncertainty']].groupby('Position').mean().reset_index()

    Ng = 20;
    block_mm_map = np.empty((Ng, Ng))
    block_mm_map[:] = np.nan
    for n in range(len(block_mm_data_mean)):
        pos = block_mm_data_mean['Position'].iloc[n]
        Positions = list(map(lambda x: int(x) - 1, pos.split('|')))
        block_mm_map[Positions[0], Positions[1]] = block_mm_data_mean['Delta ABA (kBT)'].iloc[n]

    if Plot:
        plt.figure()
        sns.heatmap(block_mm_map, cmap='coolwarm', cbar=True, vmin=0, vmax=block_mm_data_mean['Delta ABA (kBT)'].max())

        plt.xlabel('Block end', fontsize=15)
        plt.ylabel('Block start', fontsize=15)
        plt.title('$\Delta$ABA - ' + data_name, fontsize=15)
        if not Canonical:
            plt.title('$\Delta$ABA - ' + data_name + ' - Noncanonical', fontsize=15)
        ax = plt.gca()
        ax.set_xticklabels(list(map(lambda x: str(int(x)), ax.get_xticks() + 0.5)),fontsize=15);
        ax.set_yticklabels(list(map(lambda x: str(int(x)), 20 -ax.get_yticks() + 0.5)), fontsize=15,rotation=0);
        if SaveFigures:
            plt.savefig(figure_names[0],format='pdf',bbox_inches='tight')


    block_start_mm_data = block_mm_data_mean.set_index('Position').groupby(
        lambda x: int(x.split('|')[0])).mean().reset_index().rename(columns={'index': 'Block start'})

    if Plot:
        block_start_mm_data.plot(x='Block start', y='Delta ABA (kBT)', yerr='Uncertainty', marker='o', linewidth=1,
                                 color='blue')
        plt.xlabel('Block start', fontsize=15)
        plt.ylabel('$\Delta$ABA', fontsize=15)
        plt.title(data_name, fontsize=15)
        if not Canonical:
            plt.title(data_name + ' - Noncanonical', fontsize=15)


        plt.xticks(fontsize=15);
        plt.yticks(fontsize=15);
        if SaveFigures:
            plt.savefig(figure_names[1],format='pdf',bbox_inches='tight')


    return block_mm_data_mean, block_mm_map, block_start_mm_data


def plot_mut_PAM_ABA(data, data_name='Finkelstein Data', Plot=True):
    data['Mutation Type'].fillna('', inplace=True)
    data['Alignment'].fillna('', inplace=True)
    NonCanonical = data[(data['Canonical'] == False) & (data['Alignment'] != '') & (data['Mutation Type'] == '')][
        ['PAM', 'Delta ABA (kBT)', 'Uncertainty']].groupby('PAM').mean()

    if Plot:
        bar_width = 0.35
        plt.figure(figsize=(20, 6))
        plt.bar(np.arange(len(NonCanonical)), NonCanonical['Delta ABA (kBT)'], width=bar_width)
        ax = plt.gca()
        ax.set_xticks(np.arange(len(NonCanonical)) + 0.5 * bar_width);
        ax.set_xticklabels(NonCanonical.index, rotation='vertical');
        plt.xlabel('PAM', fontsize=15)
        plt.ylabel('$\Delta$ABA', fontsize=15)
        plt.title(data_name, fontsize=15)

    return NonCanonical


def predict_block_mismatches(parameters, model_id, T=60 * 10, guide_length=20, show_plot=True, show_data=True,
                             data_file = '../Data_ABA_Finkelsteinlab/champ-cas9-cas12a-data/cas9-target-e-replicate-1-delta-abas-processed.csv'):
    concentrations = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300]) #2 ** np.array(range(0, 11)) * 0.5
    reference_conc = 1 #10
    ontarget_ABA = CalcABA.calc_ABA(parameters, concentrations, reference_conc,
                                    mismatch_positions=[],
                                    model_id=model_id,
                                    guide_length=20,
                                    T=60 * 10)

    delta_ABA_mat = np.zeros((guide_length, guide_length))
    delta_ABA_mat[:] = np.nan

    for start_of_block in range(1, guide_length + 1):
        for end_of_block in range(start_of_block + 2, guide_length + 1):
            mm_block = list(range(start_of_block, end_of_block)) # Python 3! for python 2, remove the list around it!
            delta_ABA_mat[start_of_block - 1, end_of_block - 1] = CalcABA.calc_delta_ABA(parameters, concentrations,
                                                                                         reference_conc,
                                                                                         mismatch_positions=mm_block,
                                                                                         model_id=model_id,
                                                                                         guide_length=guide_length,
                                                                                         T=T,
                                                                                         ontarget_ABA=ontarget_ABA)
    if show_plot:
        ax = sns.heatmap(delta_ABA_mat, cmap='coolwarm', vmin=0, vmax=2.5)
        plt.grid()
        ax.set_xticklabels(list(map(lambda x: str(int(x)), ax.get_xticks() + 0.5)), fontsize=15);
        ax.set_yticklabels(list(map(lambda x: str(int(x)), 20 - ax.get_yticks() + 0.5)), fontsize=15, rotation=0);
        plt.xlabel('Block end', fontsize=15)
        plt.ylabel('Block start', fontsize=15)
        plt.title('Prediction', fontsize=15)

        if show_data:
            plt.figure()
            IlyaData = data_file #pd.read_csv( data_file )
            _, Blocks_heatmap, _ = plot_block_mm_ABA(IlyaData,
                                                            Plot=False,
                                                            SaveFigures=False,
                                                            figure_names=[''])
            sns.heatmap(Blocks_heatmap, cmap='coolwarm', vmin=0, vmax=2.5)
            plt.grid()
            plt.xlabel('Block end', fontsize=15)
            plt.ylabel('Block start', fontsize=15)
            plt.title('Data', fontsize=15)

    return delta_ABA_mat



def predict_1D_mmblocks(parameters, model_id, T=60 * 10, guide_length=20, show_plot=True, show_data=True,
                        data_file='../Data_ABA_Finkelsteinlab/champ-cas9-cas12a-data/cas9-target-e-replicate-1-delta-abas-processed.csv'):
    concentrations = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300]) #2 ** np.array(range(0, 11)) * 0.5
    reference_conc = 1 #10
    ontarget_ABA = CalcABA.calc_ABA(parameters, concentrations, reference_conc,
                                    mismatch_positions=[],
                                    model_id=model_id,
                                    guide_length=20,
                                    T=60 * 10)

    delta_ABA_mat = np.zeros(guide_length - 2)
    for start_of_block in range(1, guide_length - 1):
        avg_delta_ABA = 0
        count = 0
        for end_of_block in range(start_of_block + 2, guide_length + 1):
            mm_block = list(range(start_of_block, end_of_block))
            avg_delta_ABA += CalcABA.calc_delta_ABA(parameters, concentrations, reference_conc,
                                                    mismatch_positions=mm_block,
                                                    model_id=model_id,
                                                    guide_length=guide_length,
                                                    T=T,
                                                    ontarget_ABA=ontarget_ABA)
            count += 1

        delta_ABA_mat[start_of_block - 1] = avg_delta_ABA / float(count)

    if show_plot:
        ax = plt.plot(range(1, guide_length - 1),
                      delta_ABA_mat,
                      marker='o',
                      markersize=4,
                      markerfacecolor='white',
                      markeredgewidth=2,
                      linestyle='solid',
                      label='prediction model')

        sns.despine()
        plt.xticks(range(1, 19), fontsize=15);
        plt.yticks(fontsize=15);
        plt.xlabel('start of mismatched block', fontsize=15);
        plt.ylabel(r'$\Delta \rm{ABA} \ (k_BT)$', fontsize=15)

        if show_data:
            IlyaData = data_file #pd.read_csv(data_file)
            _, _, ABA_first_mm_pos = plot_block_mm_ABA(IlyaData,
                                                              Plot=False,
                                                              SaveFigures=False,
                                                              figure_names=[]);

            plt.errorbar(x=ABA_first_mm_pos['Block start'],
                         y=ABA_first_mm_pos['Delta ABA (kBT)'],
                         yerr=ABA_first_mm_pos['Uncertainty'],
                         marker='o',
                         markersize=4,
                         markerfacecolor='white',
                         markeredgewidth=2,
                         linestyle='solid',
                         capsize=5,
                         label='data Finkelstein lab'
                         )
            plt.legend(fontsize=15, loc='best')
    return delta_ABA_mat


def predict_single_mm(parameters, model_id, T=60 * 10, guide_length=20, show_plot=True, show_data=True,
                      data_file='../Data_ABA_Finkelsteinlab/champ-cas9-cas12a-data/cas9-target-e-replicate-1-delta-abas-processed.csv'):
    concentrations = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300]) #2 ** np.array(range(0, 11)) * 0.5
    reference_conc = 1 #10
    ontarget_ABA = CalcABA.calc_ABA(parameters, concentrations, reference_conc,
                                    mismatch_positions=[],
                                    model_id=model_id,
                                    guide_length=20,
                                    T=60 * 10)

    delta_ABA = np.zeros(guide_length)
    for mm_pos in range(1, guide_length + 1):
        delta_ABA[mm_pos - 1] = CalcABA.calc_delta_ABA(parameters, concentrations, reference_conc,
                                                       mismatch_positions=[mm_pos],
                                                       model_id=model_id,
                                                       guide_length=guide_length,
                                                       T=T,
                                                       ontarget_ABA=ontarget_ABA)

    if show_plot:
        ax = plt.plot(range(1, guide_length + 1),
                      delta_ABA,
                      marker='o',
                      markersize=4,
                      markerfacecolor='white',
                      markeredgewidth=2,
                      linestyle='solid',
                      label='prediction model')

        sns.despine()
        plt.xticks(range(1, 19), fontsize=15);
        plt.yticks(fontsize=15);
        plt.xlabel('mismatch position', fontsize=15);
        plt.ylabel(r'$\Delta \rm{ABA} \ (k_BT)$', fontsize=15)

        if show_data:
            IlyaData = data_file #pd.read_csv(data_file)
            single_mut_data_mean = plot_single_mut_ABA(data=IlyaData, Mut_type='r', Plot=False)

            plt.errorbar(x=single_mut_data_mean['Position'],
                         y=single_mut_data_mean['Delta ABA (kBT)'],
                         yerr=single_mut_data_mean['Uncertainty'],
                         marker='o',
                         markersize=4,
                         markerfacecolor='white',
                         markeredgewidth=2,
                         linestyle='solid',
                         label='data Finkelstein lab')
            plt.legend(fontsize=15, loc='best')
    return delta_ABA


def predict_double_mm(parameters, model_id, T=60 * 10, guide_length=20, show_plot=True, show_data=True,
                      data_file='../Data_ABA_Finkelsteinlab/champ-cas9-cas12a-data/cas9-target-e-replicate-1-delta-abas-processed.csv'):
    concentrations = np.array([0.1, 0.3, 1, 3, 10, 30, 100, 300]) #2 ** np.array(range(0, 11)) * 0.5
    reference_conc = 1 #10
    ontarget_ABA = 42#CalcABA.calc_ABA(parameters, concentrations, reference_conc,
                                   # mismatch_positions=[],
                                   # model_id=model_id,
                                   # guide_length=20,
                                    #T=60 * 10)

    delta_ABA_mat = np.zeros((guide_length, guide_length))
    for first_mm in range(1, guide_length + 1):
        for second_mm in range(1, guide_length + 1):
            delta_ABA_mat[first_mm - 1, second_mm - 1] = CalcABA.calc_ABA(parameters, concentrations,
                                                                                reference_conc,
                                                                                mismatch_positions=[first_mm,
                                                                                                    second_mm],
                                                                                model_id=model_id,
                                                                                guide_length=guide_length,
                                                                                T=T)
            # this should be calc_delta_ABA if we do not work wit rawABA dataset
            
    if show_plot:
        axHeatmap = sns.heatmap(delta_ABA_mat, cmap='coolwarm', mask=np.tril(delta_ABA_mat),vmin=2, vmax=5)
        plt.grid()
        ax = plt.gca()
        ax.set_xticklabels(list(map(lambda x: str(int(x)), ax.get_xticks() + 0.5)), fontsize=15);
        ax.set_yticklabels(list(map(lambda x: str(int(x)), 20 - ax.get_yticks() + 0.5)), fontsize=15, rotation=0);
        str_title = 'Prediction (top)'

        if show_data:
            IlyaData = data_file #pd.read_csv(data_file)
            _, double_mut_map = pltF.plot_double_mut_data(IlyaData, data_colname='ABA', Mut_type='r', Canonical=True, Ng=20, data_name='Data', Plot=False,logplot=False, SaveFigures=False, figure_name='./Figure.pdf')
            # data_colname = [delta ABA (kBT)]  if we use not rawABA dataset!
            #plot_double_mut_ABA(data=IlyaData, Mut_type='r', Plot=False)
            sns.heatmap(double_mut_map, cmap='coolwarm', ax=axHeatmap, vmin=2, vmax=5)
            str_title += ' / Data (bottom)'
        plt.title(str_title, fontsize=15)
        plt.xlabel('Mismatch 1', fontsize=15)
        plt.ylabel('Mismatch 2', fontsize=15)
    return delta_ABA_mat

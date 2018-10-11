import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('ticks')
current_colors = sns.color_palette()

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
        ax.set_yticklabels(map(lambda x: str(int(x)), ax.get_yticks() + 0.5));

    return double_mut_data_mean, double_mut_map


def make_block_ID(MM_pos_list):
    bi = MM_pos_list[0]
    bf = MM_pos_list[-1]
    block = np.arange(bi, bf + 1)
    if (len(block) != len(MM_pos_list)) or np.any(block != np.array(MM_pos_list)):
        return ''
    return '|'.join(map(str, [bi, bf]))


def plot_block_mm_ABA(data, data_name='Finkelstein Data', Canonical=True, Plot=True,
                      SaveFigures=False, figure_names=['figure1','figure 2']):
    data['Mutation Type'].fillna('', inplace=True)
    select_multi_mm = data['Mutation Type'].apply(lambda x: np.unique(np.array(x.split('|')))[0] == 'r')
    multi_mm_data = data[(data['Canonical'] == Canonical) & (select_multi_mm) & (data['Mutation Count'] > 2)][['Mutation ID', 'Delta ABA (kBT)', 'Uncertainty']]
    multi_mm_data['MM Positions'] = multi_mm_data['Mutation ID'].apply(
        lambda x: map(lambda y: int(y.split(':')[1]), x.split('|')))
    multi_mm_data['Position'] = multi_mm_data['MM Positions'].apply(make_block_ID)
    multi_mm_data = multi_mm_data[(multi_mm_data['Position'] != '')]
    block_mm_data_mean = multi_mm_data[['Position', 'Delta ABA (kBT)', 'Uncertainty']].groupby('Position').mean().reset_index()

    Ng = 20;
    block_mm_map = np.empty((Ng, Ng))
    block_mm_map[:] = np.nan
    for n in range(len(block_mm_data_mean)):
        pos = block_mm_data_mean['Position'].iloc[n]
        Positions = map(lambda x: int(x) - 1, pos.split('|'))
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
        ax.set_xticklabels(map(lambda x: str(int(x)), ax.get_xticks() + 0.5),fontsize=15);
        ax.set_yticklabels(map(lambda x: str(int(x)), 20 -ax.get_yticks() + 0.5), fontsize=15,rotation=0);
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
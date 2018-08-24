import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('ticks');
current_colors = sns.color_palette()
import sys
sys.path.append('../code_Boyle/')
import CRISPR_dCas9_binding_curve_Boyle as dCas9
reload(dCas9);


def load_simm_anneal(filename, Nparams):
    '''
    Load the parameter set from simmulated annealing.
    Fix indexing based on parameterisation to correctly import the table
    :param filename:
    :param Nparams:
    :return:
    '''
    fit = pd.read_csv(filename, delimiter='\t', index_col=Nparams+2)
    fit = fit.reset_index()
    final_result = []
    for param in range(1, Nparams + 1):
        col = 'Parameter ' + str(param)
        final_result.append(fit[col].iloc[-1])

    sa_result = np.array(final_result)
    return sa_result

def calc_predictions(parameters,model_id):
    '''
    Predict occupancy (12 hrs), on-rate and off-rate for all combinations of two mismachtes
    Uses 'CRISPR_dCas9_binding_curve_Boyle' module

    :param parameters: (fitted) set of relevant/free epsilon and rate values
    :param model_id: naming used to unpack parameter set/ interpret energy landscape
    :return: three matrices. Diagonals contain single-mismatches
    '''
    on_rate_predict = np.zeros([20, 20])
    off_rate_predict = np.zeros([20, 20])
    Pbound_predict = np.zeros([20, 20])
    Pbound_OT, _, _ = dCas9.calc_Boyle(True, True, True, parameters, [],model_id=model_id)
    for i in range(0, 20):
        for j in range(0, 20):
            mismatch_positions = [20 - i, 20 - j]
            if i == j:
                mismatch_positions = [20 - i]
            predictions = dCas9.calc_Boyle(True, True, True, parameters, mismatch_positions,model_id=model_id)
            Pbound_predict[i, j] = predictions[0] / Pbound_OT
            on_rate_predict[i, j] = predictions[1]
            off_rate_predict[i, j] = predictions[2]
    return Pbound_predict, on_rate_predict, off_rate_predict


def plot_heatmap(model ,kind='Occupancy', fldr_Boyle_data = '../Data_Boyle/KoenDataForMisha/BoyleData/'):
    '''
    Plot heatmap for double mismatches with Boyle's data in below diagonal and model above the diagonal
    :param model: matrix with double-mismatch values
    :param kind: 'Occupancy', 'OnRate' or 'OffRate'
    :param fldr_Boyle_data:
    :return:
    '''
    # 1) settings based on physical quantity you want to plot
    if kind == 'Occupancy':
        colormap = 'Greens'
        val_min = 0
        val_max = 1
        data_file = 'DataOccupancy.txt'
        title = 'Occupancy (after 12 hrs) \n relative to cognate site'
    elif kind == 'OnRate':
        colormap = 'Reds'
        val_min = 0
        val_max = 0.0003
        data_file = 'DataOnRate.txt'
        title = r'Associaton rate 1nM ($s^{-1}$)'
    elif kind == 'OffRate':
        colormap = 'Blues'
        val_min = 0
        val_max = 0.0002
        data_file = 'DataOffRate.txt'
        title = r'Dissociaton rate 10nM ($s^{-1}$)'

    # 2) load the experimental dataset
    experiment = np.loadtxt(delimiter=',' , fname=fldr_Boyle_data + data_file)

    # 3) plot the data:
    mask_exp = np.zeros(shape=experiment.shape)
    for i in range(len(experiment)):
        for j in range(i - 1, len(experiment)):
            mask_exp[i, j] = 1
    sns.heatmap(experiment, cmap=colormap, mask=mask_exp, cbar=True, vmin=0, vmax=None);

    # 4) Plot the model prediction alongside
    mask = np.ones(shape=model.shape)
    for i in range(len(model)):
        for j in range(i + 1, len(model)):
            mask[i, j] = 0
    sns.heatmap(model, cmap=colormap, mask=mask, cbar=True, vmin=val_min, vmax=val_max);

    # 5) Adjust ticks and labels to get correct nucleotide positions
    plt.xticks([i + 0.5 for i in range(20)],
               [20, '', '', '', 15, '', '', '', '', 10, '', '', '', '', 5, '', '', '', '', 1], rotation=0,
               fontsize=15);
    plt.yticks([i + 0.5 for i in range(20)],
               [1, '', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0,
               fontsize=15);

    # 6) Further window dressing
    plt.xlabel('mismatch 1', fontsize=15)
    plt.ylabel('mismatch 2', fontsize=15)
    plt.title(title, fontsize=15);
    return

def plot_single_mismatches(model ,kind='Occupancy', fldr_Boyle_data = '../Data_Boyle/KoenDataForMisha/BoyleData/'):
    # 1) settings based on physical quantity you want to plot
    if kind == 'Occupancy':
        color = 'green'
        val_max = 1.2
        val_min = -(0.03) * val_max
        data_file = 'DataOccupancy.txt'
        ylabel = 'Occupancy (after 12 hrs) \n relative to cognate site'
    elif kind == 'OnRate':
        color = 'red'
        val_max = 0.0003
        val_min = -(0.03)* val_max
        data_file = 'DataOnRate.txt'
        ylabel = r'Associaton rate 1nM ($s^{-1}$)'
    elif kind == 'OffRate':
        color = 'blue'
        val_max = 0.00020
        val_min = -(0.03) * val_max
        data_file = 'DataOffRate.txt'
        ylabel = r'Dissociation rate 10nM ($s^{-1}$)'

    # 2) load the experimental dataset
    experiment = np.diag( np.loadtxt(delimiter=',', fname=fldr_Boyle_data + data_file) )

    # 3) extract single-mismatches from entire prediction
    model= np.diag(model)

    #4) plot
    positions = 20 - np.arange(0,20)
    plt.plot(positions, experiment, linestyle='dashed', marker='s', color=color,label='Data')
    plt.plot(positions, model, linestyle='solid', marker='^', color=color, label='Model')

    #5) aesthetic and window dressing
    plt.ylim(val_min,val_max )
    plt.xlim(1, 21)
    plt.xlabel('mismatch position', fontsize=15)
    plt.ylabel(ylabel, fontsize=15)
    plt.xticks([i + 0.5 for i in range(20)],
               [1, '', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0,
               fontsize=15);
    plt.yticks(fontsize=15)
    plt.legend(loc='best', fontsize=12, frameon=True)
    sns.despine()
    return

def plot_landscape(parameters, model_id):
    '''
    Plot the free-energy landscape of the on-target
    :param parameters:
    :param model_id:
    :return:
    '''
    epsilon, fwrd_rates = dCas9.unpack_parameters(parameters, model_id, guide_length=20)
    epsilon_C = epsilon[:21]
    epsilon_C[1:] *= -1
    landscape = [0.0]
    for eps in epsilon_C:
        landscape.append(landscape[-1] + eps)
    plt.plot(landscape, marker='s')

    # window dressing:
    plt.xlabel('targeting progression', fontsize=15)
    plt.ylabel(r'free-energy ($k_BT$)',fontsize=15)
    plt.xticks([i for i in range(1,21)],
               [1, '', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0
               ,fontsize=15);
    plt.yticks(fontsize=15)
    sns.despine()
    return landscape


def plot_mismatch_penalties(parameters, model_id):
    '''
    plot mismatch penalties VS position as a bar plot
    :param parameters:
    :param model_id:
    :return:
    '''
    epsilon, fwrd_rates = dCas9.unpack_parameters(parameters, model_id, guide_length=20)
    epsilon_I = epsilon[21:]
    plt.bar([i + 1 for i in range(20)], epsilon_I)

    # window dressing:
    plt.xlabel('targeting progression', fontsize=15)
    plt.ylabel(r'free-energy ($k_BT$)',fontsize=15)
    plt.xticks([i for i in range(1,21)],
               [1, '', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0
               ,fontsize=15);
    plt.yticks(fontsize=15)
    return
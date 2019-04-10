import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('ticks');
current_colors = sns.color_palette()
import sys
import importlib as imp
sys.path.append('../code_Boyle/')
import CRISPR_dCas9_binding_curve_Boyle as dCas9
imp.reload(dCas9)
sys.path.append('../code_general/')
import read_model_ID
imp.reload(read_model_ID);
import Calculate_ABA_Finkelsteinlab_Diewertje as CalcABA

def load_simm_anneal(filename, Nparams, fatch_solution='final'):
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

        if fatch_solution == 'final':
            final_result.append(fit[col].iloc[-1])
        else:
            final_result.append(fit[col].iloc[fatch_solution])

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
    Pbound_OT, _, _ = calc_Boyle(True, True, True, parameters, [],model_id=model_id)
    for i in range(0, 20):
        for j in range(0, 20):
            mismatch_positions = [20 - i, 20 - j]
            if i == j:
                mismatch_positions = [20 - i]
            predictions = calc_Boyle(True, True, True, parameters, mismatch_positions,model_id=model_id)
            Pbound_predict[i, j] = predictions[0] / Pbound_OT
            on_rate_predict[i, j] = predictions[1]
            off_rate_predict[i, j] = predictions[2]
    return Pbound_predict, on_rate_predict, off_rate_predict


def plot_heatmap(model ,kind='Occupancy', fldr_Boyle_data = '../Data_Boyle/KoenDataForMisha/BoyleData/',show_plot=True,axis=None,cbar=True):
    '''
    Plot heatmap for double mismatches with Boyle's data in below diagonal and model above the diagonal
    :param model: matrix with double-mismatch values
    :param kind: 'Occupancy', 'OnRate' or 'OffRate'
    :param fldr_Boyle_data:
    :return:
    '''

    Ng = 20
    # 1) settings based on physical quantity you want to plot
    if kind == 'Occupancy':
        colormap = 'Greens'
        data_file = 'DataOccupancy.txt'
        title = 'Occupancy (after 12 hrs) \n relative to cognate site'
    elif kind == 'OnRate':
        colormap = 'Reds'
        data_file = 'DataOnRate.txt'
        title = r'Associaton rate 1nM ($s^{-1}$)'
    elif kind == 'OffRate':
        colormap = 'Blues'
        data_file = 'DataOffRate.txt'
        title = r'Dissociaton rate 10nM ($s^{-1}$)'


    # 2) load the experimental dataset
    experiment = np.loadtxt(delimiter=',' , fname=fldr_Boyle_data + data_file)

    val_min = 0
    val_max_1 = np.nanmax(experiment)
    val_max_2 = np.max(model)
    val_max = np.max(np.array([val_max_1, val_max_2]))
    val_max = (1+0.05)*val_max
    if kind == 'Occupancy':
        val_max = 1

    if show_plot:
        if axis is None:
            ax = plt.gca()
        else:
            ax = axis

    # 3) plot the data:
    mask_exp = np.zeros(shape=experiment.shape)
    for i in range(len(experiment)):
        for j in range(i - 1, len(experiment)):
            mask_exp[i, j] = 1
    if show_plot:
        sns.heatmap(experiment, cmap=colormap, mask=mask_exp, cbar=cbar, vmin=0, vmax=val_max,ax=ax);

    # 4) Plot the model prediction alongside
    mask = np.ones(shape=model.shape)
    for i in range(len(model)):
        for j in range(i + 1, len(model)):
            mask[i, j] = 0
    if show_plot:
        sns.heatmap(model, cmap=colormap, mask=mask, cbar=cbar, vmin=val_min, vmax=val_max,ax=ax);


        # 5) Adjust ticks and labels to get correct nucleotide positions
        ax.set_xticklabels(map(lambda x: str(int(Ng-x)), ax.get_xticks() - 0.5));
        #ax.set_yticklabels(map(lambda x: str(int(Ng-x)), ax.get_yticks() - 0.5));
        ax.set_yticklabels(map(lambda x: str(int(x+1)), ax.get_yticks() - 0.5));

        # 6) Further window dressing
        ax.set_xlabel('mismatch 1', fontsize=15)
        ax.set_ylabel('mismatch 2', fontsize=15)
        ax.set_title(title, fontsize=15);
    return model, experiment

def plot_single_mismatches(model ,kind='Occupancy', fldr_Boyle_data = '../Data_Boyle/KoenDataForMisha/BoyleData/',show_plot=True, axis=None):
    # 1) settings based on physical quantity you want to plot
    if kind == 'Occupancy':
        color = 'green'
        data_file = 'DataOccupancy.txt'
        ylabel = 'Occupancy (after 12 hrs) \n relative to cognate site'
    elif kind == 'OnRate':
        color = 'red'
        data_file = 'DataOnRate.txt'
        ylabel = r'Associaton rate 1nM ($s^{-1}$)'
    elif kind == 'OffRate':
        color = 'blue'
        data_file = 'DataOffRate.txt'
        ylabel = r'Dissociation rate 10nM ($s^{-1}$)'

    # 2) load the experimental dataset
    experiment = np.diag( np.loadtxt(delimiter=',', fname=fldr_Boyle_data + data_file) )

    # 3) extract single-mismatches from entire prediction
    model= np.diag(model)

    val_max_1 = np.nanmax(experiment)
    val_max_2 = np.max(model)
    val_max = np.max(np.array([val_max_1, val_max_2]))
    val_max = (1+0.05)*val_max
    val_min = -(0.03) * val_max

    #4) plot
    if show_plot:
        positions = 20 - np.arange(0,20)
        if axis is None:
            plt.plot(positions, experiment, linestyle='dashed', marker='s', color=color,label='Data')
            plt.plot(positions, model, linestyle='solid', marker='^', color=color, label='Model')
            ax=plt.gca()
        else:
            ax = axis
            ax.plot(positions, experiment, linestyle='dashed', marker='s', color=color, label='Data')
            ax.plot(positions, model, linestyle='solid', marker='^', color=color, label='Model')


        #5) aesthetic and window dressing
        ax.set_ylim(val_min,val_max )
        ax.set_xlim(1, 21)
        ax.set_xlabel('mismatch position', fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)
        ax.set_xticks([i + 0.5 for i in range(20)])
        ax.set_xticklabels([1, '', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20],
                           rotation=0,
                           fontsize=15);
        ax.set_yticklabels(ax.get_yticks(),fontsize=15)
        ax.legend(loc='best', fontsize=12, frameon=True)
        sns.despine(ax=ax)
    return model, experiment


##### Extra function for calculate predicitons from the file: 
#### CRISPR_dCas9_binding_curve_Boyle 
def calc_Boyle(CalcOccupancy, CalcOffRate, CalcOnRate,
               parameters, mismatch_positions, guide_length=20, model_id='general_energies'):
    '''
    :param CalcOccupancy: boolian True if you want to calculate occupancy after 12 hrs
    :param CalcOffRate: boolian, True if you want to calculate dissocation rate
    :param CalcOnRate:  boolian, True if you want to calculate association rate
    :return:
    '''
    # 0) Initialize output (in case I do not want to calculate all 3:
    bound_fraction = float('nan')
    dissociation_rate = float('nan')
    association_rate = float('nan')

    # 1) Unpack parameters and build rate matrix to be used in all types of calculations
    rate_matrix = CalcABA.get_master_equation(parameters, mismatch_positions, model_id, guide_length)

    #2) Do I want (need) the occupancy after 12 hours?
    if CalcOccupancy or CalcOffRate:
        everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))
        Probability = CalcABA.get_Probability(rate_matrix=rate_matrix,initial_condition=everything_unbound,T=12*3600)
        bound_fraction = np.sum(Probability[1:])

    #3) Do I want the dissociation rate ?
    if CalcOffRate:
        dissociation_rate = calc_dissociation_rate(rate_matrix=rate_matrix,
                                                 initial_condition=Probability,timepoints=[500.,1000.,1500.])
    #4) Do I want the association rate?
    if CalcOnRate:
        association_rate = calc_association_rate(rate_matrix=rate_matrix,timepoints=[500.,1000.,1500.],
                                                 guide_length=guide_length)

    return bound_fraction, association_rate, dissociation_rate

def calc_association_rate(rate_matrix,timepoints=[500.,1000.,1500.],guide_length=20,rel_concentration=0.1):
    '''
    Association experiment for the effective on-rate:
    (repeating the protocol as used in Boyle et al.)
    :return:
    '''

    #1) Calculate bound fraction for specified time points:
    # Association experiment starts with all dCas9 being unbound:
    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))

    #3) Association rate is taken at 1nM the other data at 10nM --> adjust k_on appropriatly:
    new_rate_matrix = rate_matrix.copy()
    new_rate_matrix[0][0] *= rel_concentration  #rel_concentration=1 corresponds to 10nM
    new_rate_matrix[1][0] *= rel_concentration

    bound_fraction = []
    for time in timepoints:
        Probabilities = CalcABA.get_Probability(rate_matrix=new_rate_matrix,initial_condition=everything_unbound,T=time)
        bound_fraction.append( np.sum(Probabilities[1:]))

    #2) Use least-squares to fit straight line with origin forced through zero:
    kon_lin_fit = least_squares_line_through_origin(x_points=np.array(timepoints),y_points=np.array(bound_fraction))
    return kon_lin_fit

def calc_dissociation_rate(rate_matrix,initial_condition,
                           timepoints=[500.,1000.,1500.]):
    '''
    Dissociation experiment for (effective off-rate):
    method as used in Boyle et al.

    within calc_Boyle() (main function call) you will all ready calculate the occupation before flushing out dCas9
    :return:
    '''

    # 1) Flush out all freely floating dCas9 --> Set occupation in solution to zero:
    after_12hrs = initial_condition.copy()
    after_12hrs[0] = 0.0

    #2) Renormalize remainder:
    after_12hrs = after_12hrs / np.sum(after_12hrs)

    #3) Flush out all freely floating dCas9 --> Set on-rate to zero and rebuild matrix:
    new_rate_matrix = rate_matrix.copy()
    new_rate_matrix[0][0] = 0.0
    new_rate_matrix[1][0] = 0.0

    #4)  For time_k in timepoints solve the Master equation and track the fraction of dCas9 molecules in solution:
    unbound_fraction = []
    for time in timepoints:
        Probabilities = CalcABA.get_Probability(rate_matrix=new_rate_matrix, initial_condition=after_12hrs, T=time)
        unbound_fraction.append(Probabilities[0])

    # 5) Use least-squares to fit :
    _ , k_off_fit = least_squares(x_points=timepoints, y_points=unbound_fraction)
    return k_off_fit

def least_squares(x_points, y_points):
    size = len(x_points)
    X = np.ones((size, 2))
    X[:, 1] = x_points

    XT = X.transpose()
    Y = y_points

    a = np.dot(np.dot(np.linalg.inv(np.dot(XT, X)), XT), Y)

    intercept = a[0]
    slope = a[1]
    return intercept, slope

def least_squares_line_through_origin(x_points, y_points):
    return np.sum( x_points*y_points )/np.sum( x_points*x_points )

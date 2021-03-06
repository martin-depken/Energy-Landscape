import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style('ticks')
import read_model_ID
import importlib as imp
imp.reload(read_model_ID);



def plot_free_energy_landscape(parameters,model_id,show_plot=True):
    '''
    Plot the approximate the free-energy landscape for the on-target

    Note: Technically, the system just has a single free-energy:
        F = -log( sum_i e^{-E_i})
    We will get F_n by  combining the Free-Energies for landscapes with mismatches from (n+1) on
    In this way, we may assume states after the block of mismatches to be inaccessible (the regime the data dictates)
    Hence, we can get an approximate free-energy landscape for the on-target in this manner.

    :param parameters:
    :param model_id:
    :return:
    '''
    epsilon, fwrd_rates = read_model_ID.unpack_parameters(parameters, model_id, guide_length=20)

    Epsilon = epsilon.copy()

    epsilon_C = np.array(Epsilon[:21])
    epsilon_C[1:] *= -1
    epsilon_C=np.ndarray.tolist(epsilon_C)

    landscape = [0.0]
    for eps in epsilon_C:
        landscape.append(landscape[-1] + eps)
    landscape = np.array(landscape)

    if show_plot:
        plt.figure()
        plt.plot(range(-1,21),landscape, marker='s');
        plt.xlabel('targeting progression', fontsize=15)
        plt.ylabel(r'(approx.) free-energy ($k_BT$)',fontsize=15)
        plt.xticks(range(-1,21),
                   [ 'S','P',1,'', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0
                   ,fontsize=15);
        plt.yticks(fontsize=15)
        plt.grid('on')
        sns.despine()




    free_energy =  -1*np.log( np.cumsum( np.exp(-landscape[1:])   )  )

    # plotting:
    if show_plot:
        plt.figure()
        plt.plot(range(0,21),free_energy, marker='s');
        plt.xlabel('targeting progression', fontsize=15)
        plt.ylabel(r'(approx.) free-energy ($k_BT$)',fontsize=15)
        plt.xticks(range(0,21),
                   [ 'P',1,'', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0
                   ,fontsize=15);
        plt.yticks(fontsize=15)
        plt.grid('on')
        sns.despine()

    return landscape, free_energy



def plot_landscape(parameters, model_id, axis=None):
    '''
    Plot the (free-)energy landscape of the on-target
    :param parameters:
    :param model_id:
    :return:
    '''
    epsilon, fwrd_rates = read_model_ID.unpack_parameters(parameters, model_id, guide_length=20)
    epsilon_C = epsilon[:21]
    epsilon_C[1:] *= -1

    landscape = [0.0]
    for eps in epsilon_C:
        landscape.append(landscape[-1] + eps)

    if axis is None:
        plt.plot(range(-1,21),landscape, marker='s');
        ax = plt.gca()
    else:
        ax = axis
        ax.plot(range(-1,21),landscape, marker='s');


    # window dressing:
    ax.set_xlabel('targeting progression', fontsize=15)
    ax.set_ylabel(r'free-energy ($k_BT$)',fontsize=15)
    ax.set_xticks(range(-1,21))
    ax.set_xticklabels( ['S', 'P',1,'', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20],
                        rotation=0,
                        fontsize=15);
    ax.set_yticklabels(ax.get_yticks() ,fontsize=15)
    ax.grid('on')
    sns.despine(ax=ax)
    return landscape


def plot_mismatch_penalties(parameters, model_id,axis=None):
    '''
    plot mismatch penalties VS position as a bar plot
    :param parameters:
    :param model_id:
    :return:
    '''
    epsilon, fwrd_rates = read_model_ID.unpack_parameters(parameters, model_id, guide_length=20)
    epsilon_I = epsilon[21:]

    if axis is None:
        ax = plt.gca()
    else:
        ax = axis
    ax.bar(range(1,21), epsilon_I, alpha=0.5)
    # window dressing:
    ax.set_xlabel('targeting progression', fontsize=15)
    ax.set_ylabel(r'mismatch penalties ($k_BT$)',fontsize=15)
    ax.set_xticks(np.arange(1,21)+0.5)
    ax.set_xticklabels([1, '', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20],
                       rotation=0,
                       fontsize=15);
    ax.set_yticklabels(ax.get_yticks(),fontsize=15)
    return epsilon_I
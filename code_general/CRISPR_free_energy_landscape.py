import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
sns.set_style('ticks')
import read_model_ID
reload(read_model_ID);



def plot_landscape(parameters, model_id):
    '''
    Plot the free-energy landscape of the on-target
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
    plt.plot(range(-1,21),landscape, marker='s')

    # window dressing:
    plt.xlabel('targeting progression', fontsize=15)
    plt.ylabel(r'free-energy ($k_BT$)',fontsize=15)
    plt.xticks(range(-1,21),
               ['S', 'P',1,'', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0
               ,fontsize=15);
    plt.yticks(fontsize=15)
    plt.grid('on')
    sns.despine()
    return landscape


def plot_mismatch_penalties(parameters, model_id):
    '''
    plot mismatch penalties VS position as a bar plot
    :param parameters:
    :param model_id:
    :return:
    '''
    epsilon, fwrd_rates = read_model_ID.unpack_parameters(parameters, model_id, guide_length=20)
    epsilon_I = epsilon[21:]
    plt.bar(range(1,21), epsilon_I)

    # window dressing:
    plt.xlabel('targeting progression', fontsize=15)
    plt.ylabel(r'mismatch penalties ($k_BT$)',fontsize=15)
    plt.xticks(np.arange(1,21)+0.5,
               [1, '', '', '', 5, '', '', '', '', 10, '', '', '', '', 15, '', '', '', '', 20], rotation=0
               ,fontsize=15);
    plt.yticks(fontsize=15)
    return epsilon_I
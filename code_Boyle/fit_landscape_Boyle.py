import numpy as np
import functools

import Boyle_data_processing
import CRISPR_dCas9_binding_curve_Boyle as CRISPR
import SimulatedAnnealing_Boyle_parallel as SA


'''
********************
Main script to be run on cluster for Simmulated Annealing fitting of energy landscape to data from Boyle et al.

Misha Klein and Behrouz Eslami-Mossallam
********************

1. Load the datafiles from Boyle and preprocess it such that it allows for multiprocessing 
2. feed data and model into SA code 
'''
#################
# /* Settings *\#
#################

# Loading the data
path_to_Boyle_data='../Data_Boyle/'
replica_ID = '1'


# Simmulated Annealing
nmbr_cores = 1
model_ID = 'general_energies'
monitor_file = 'monitor.txt'
fit_result_file = 'fit.txt'


gRNA_length = 20
Weights_Datasets = Boyle_data_processing.weights_averages(replica_ID,path_to_Boyle_data)
upper_bnd = [100.0]*41 +[100.0] + [1000.0]  # estimated upper bounds based on Koen's previous results.
lower_bnd = [0.0]*43    #last element is rate from solution to PAM. Second to last is internal forward rate
initial_guess = [5.0]*41 + [1.0] + [100.]


###########################
# /* Objective function *\#
###########################
KineticModel = functools.partial(CRISPR.calc_Chi_square,
                    weights=Weights_Datasets,
                    guide_length=gRNA_length,
                    model_id=model_ID)



##############################################
# /* Preprocess the data from Boyle et al. *\#
##############################################
xdata, ydata = Boyle_data_processing.prepare_multiprocessing(replica_ID,path_to_Boyle_data)





##############################################
# /*   Call the Simulated Annealing code   *\#
##############################################
fit_result = SA.sim_anneal_fit(xdata=xdata,
                               ydata=ydata,
                               yerr = [],
                               Xstart= np.array(initial_guess),
                               lwrbnd= np.array(lower_bnd),
                               upbnd= np.array(upper_bnd),
                            model='I_am_using_multi_processing_in_stead',
                            objective_function=KineticModel,
                            Tstart=0.1,
                            delta=np.log(2.0),
                            tol=1E-3,
                            Tfinal=0.0,
                            adjust_factor=1.1,
                            cooling_rate=0.95,
                            N_int=1000,
                            AR_low=40,
                            AR_high=60,
                            use_multiprocessing=True,
                            nprocs=nmbr_cores,
                            use_relative_steps=True,
                            output_file_results = fit_result_file,
                            output_file_monitor = monitor_file
                               )









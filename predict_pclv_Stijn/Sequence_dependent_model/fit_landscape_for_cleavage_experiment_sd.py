import numpy as np
import functools
import sys
import multiprocessing as mp

PATH_HPC05 = '/home/svandersmagt/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)
import preprocessing as processing
import calculate_cleavage_rate_sd as CRISPR
import SimulatedAnnealing_Nucleaseq_parallel_sd as SA

from time import time
'''
********************
Main script to be run on cluster for Simmulated Annealing fitting of energy landscape to data from Boyle et al.

Misha Klein and Behrouz Eslami-Mossallam
********************

1. Load the datafiles from Boyle and preprocess it such that it allows for multiprocessing 
2. feed data and model into SA code 
'''


def main(argv):
    #################
    # /* Settings *\#
    #################
    # Loading the data
    use_cluster = bool(int(argv[8]))
    if use_cluster:
        path_to_dataClv= PATH_HPC05+'data_nucleaseq_Finkelsteinlab/targetE/'
        nmbr_cores = 19
    else:
        # On my laptop use this line in stead:
        path_to_dataClv = '../../' + '/data_nucleaseq_Finkelsteinlab/targetE/'
        nmbr_cores = 2

    # Simmulated Annealing
    combined_fit = bool(int(argv[9]))
    
    filename_clv = argv[1]
    model_ID = argv[3]
    monitor_file = argv[5]
    fit_result_file = argv[6]
    init_monitor_file = argv[7]
    fit_to_wa = True
    
    gRNA_length = 20
    #fit_to_median = False   
    
    upper_bnd =      [4.]*33  + [10.]*20
    lower_bnd =      [-4.]*33 + [0.]*20
    initial_guess =  [0.]*33  + [5.]*20
    

    ###########################
    # /* Objective function *\#
    ###########################
    
    KineticModel = functools.partial(CRISPR.calc_chi_squared,
                        guide_length=gRNA_length,
                        model_id=model_ID)


    #############################################
    # /* Preprocess the data from Boyle et al. *\#
    ##############################################
    if fit_to_wa:
        xdata, ydata, yerr = processing.prepare_multiprocessing_seq_dep_wa(filename_clv,path_to_dataClv)
    else:
        xdata, ydata, yerr = processing.prepare_multiprocessing_seq_dep(filename_clv,path_to_dataClv)
    
     
    ##############################################
    # /*   Call the Simulated Annealing code   *\#
    ##############################################



    t1 = time()
    fit_result = SA.sim_anneal_fit(xdata=xdata,
                                   ydata=ydata,
                                   yerr = yerr,
                                   Xstart= np.array(initial_guess),
                                   lwrbnd= np.array(lower_bnd),
                                   upbnd= np.array(upper_bnd),
                                model='I_am_using_multi_processing_in_stead',
                                objective_function=KineticModel,
                                Tstart=10000.,             # infered from run on my computer/other runs on cluster
                                use_relative_steps=False,
                                delta=0.1,
                                tol=1E-5,
                                Tfinal=0.0,
                                potential_threshold = 500.,
                                adjust_factor=1.1,
                                cooling_rate_high=0.95,
                                cooling_rate_low=0.99,
                                N_int=1000,
                                Ttransition=1000.,
                                AR_low=40,
                                AR_high=60,
                                use_multiprocessing=True,
                                nprocs=nmbr_cores,
                                output_file_results = fit_result_file,
                                output_file_monitor = monitor_file,
                                output_file_init_monitor=init_monitor_file,
                                NMAC=False, #non-monotonic adaptive cooling
                                reanneal=False, #reheating when in local minimum, set to False to do no reheating
                                random_step=False
                                )

    t2 = time()

    print "Time elapsed: " + str((t2-t1)/(3600.)) + ' hrs'

    return






if __name__ == "__main__":
    main(sys.argv)


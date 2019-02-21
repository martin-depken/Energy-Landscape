import numpy as np
import functools
import sys

PATH_HPC05 = '/home/svandersmagt/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)
import Nucleaseq_data_processing as processing
import calculate_cleavage_rate as CRISPR
import SimulatedAnnealing_Nucleaseq_parallel as SA

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
    use_cluster = bool(int(argv[6]))
    if use_cluster:
        path_to_data= PATH_HPC05+'/data_nucleaseq_Finkelsteinlab/targetD/'
        nmbr_cores = 19
    else:
        # On my laptop use this line in stead:
        path_to_data = '../' + '/data_nucleaseq_Finkelsteinlab/targetD/'
        nmbr_cores = 1

    # Simmulated Annealing
    filename = argv[1]
    model_ID =  argv[2]
    monitor_file = argv[3]
    fit_result_file = argv[4]
    init_monitor_file = argv[5]
    
    gRNA_length = 20
    times = [0,12,60,180,600,1800,6000,18000,60000]
##### This should change according to the model.
    upper_bnd =  [10.0] + [10.0]*40 +  [3.0] *2
    lower_bnd = [0.0] + [-10.0]*40 + [-7.0] *2
    initial_guess =  [5.0] + [0.0]*40 + [0.0] *2

    ###########################
    # /* Objective function *\#
    ###########################
    KineticModel = functools.partial(CRISPR.calc_chi_squared,
                        times = times,
                        guide_length=gRNA_length,
                        model_id=model_ID)


    #############################################
    # /* Preprocess the data from Boyle et al. *\#
    ##############################################
    xdata, ydata, yerr = processing.prepare_multiprocessing_nucleaseq(filename,path_to_data)
    # print ydata

    # print "test ... " + ' \n'
    # KineticModel(np.array(initial_guess),xdata,ydata,np.array([]),1.0)


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
                                Tstart=100.,             # infered from run on my computer/other runs on cluster
                                use_relative_steps=False,
                                delta=1.0,
                                tol=1E-5,
                                Tfinal=0.0,
                                adjust_factor=1.1,
                                cooling_rate=0.99,
                                N_int=1000,
                                AR_low=40,
                                AR_high=60,
                                use_multiprocessing=True,
                                nprocs=nmbr_cores,
                                output_file_results = fit_result_file,
                                output_file_monitor = monitor_file,
                                output_file_init_monitor=init_monitor_file
                                   )

    t2 = time()

    print "Time elapsed: " + str((t2-t1)/(3600.)) + ' hrs'

    return






if __name__ == "__main__":
    main(sys.argv)


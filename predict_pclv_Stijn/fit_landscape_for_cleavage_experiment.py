import numpy as np
import functools
import sys
import multiprocessing as mp

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
    use_cluster = bool(int(argv[8]))
    if use_cluster:
        path_to_dataClv= PATH_HPC05+'data_nucleaseq_Finkelsteinlab/targetE/'
        path_to_dataAba= path_to_dataClv
        nmbr_cores = 19
    else:
        # On my laptop use this line in stead:
        path_to_dataClv = '../' + '/data_nucleaseq_Finkelsteinlab/targetE/'
        path_to_dataAba = '../data_ABA_Finkelsteinlab/champ-cas9-cas12a-data/'
        nmbr_cores = 2

    # Simmulated Annealing
    combined_fit = bool(int(argv[9]))
    
    filename_clv = argv[1]
    filename_aba = argv[2]
    model_ID_aba = argv[4]
    if combined_fit:
        model_ID =  argv[3] + '+' + argv[4]
    else:
        model_ID = argv[3]
    monitor_file = argv[5]
    fit_result_file = argv[6]
    init_monitor_file = argv[7]
    
    #----------
    # FILL IN:
    fit_to_wa = True #Fit to the weighted average instead of full data
    only_single = False #Fit only to single mismatches and on-target data
    combined_boyle = False #Combine Nucleaseq data with Boyle data
    log_on_boyle = False #Use logarithmic values when fitting to Boyle data
    fill_in_boyle = False #Fill in missing values with estimates when fitting to Boyle data
    combined_CHAMP = True #Combine Nucleaseq data with CHAMP data
    
    gRNA_length = 20
    
    upper_bnd =      [5.]  + [10.]*20   + [10.]*20 + [5.]    + [3.]  + [4.]
    lower_bnd =      [-1.]  + [-10.]*20  + [0.]*20  + [-5.]   + [1.]  + [-1.]
    initial_guess =  [2.]  + [0.]*20    + [5.]*20  + [0.]  + [2.]  + [2.]
    #----------

    ###########################
    # /* Objective function *\#
    ###########################
    
    #concentrations = np.array([1.,30.,100.])
    concentrations = np.array([0.1, 0.3, 1., 3., 10., 30., 100., 300.]) # in nanoMolair
    reference = 10. # in nanomolair, important: use float, not int
    
    KineticModel = functools.partial(CRISPR.calc_chi_squared,
                        guide_length=gRNA_length,
                        model_id=model_ID,
                        log_on=log_on_boyle,combined_boyle=combined_boyle, combined_CHAMP=combined_CHAMP,
                        concentrations=concentrations, reference=reference)


    #############################################
    # /* Preprocess the data from Boyle et al. *\#
    ##############################################
    if combined_fit:
        if combined_CHAMP:
            xdata, ydata, yerr = processing.prepare_multiprocessing_combined_aba(filename_aba,filename_clv,path_to_dataAba,path_to_dataClv,fit_to_wa)
        if combined_boyle:
            xdata, ydata, yerr = processing.prepare_multiprocessing_combined('1',filename_clv,path_to_dataAba,path_to_dataClv,fit_to_wa,only_single,log_on_boyle,fill_in_boyle)
    
    if not combined_fit:
        xdata, ydata, yerr = processing.prepare_multiprocessing_nucleaseq_log(filename_clv,path_to_dataClv,fit_to_wa,only_single)
    
    
    #######################
    # Determining weights #
    #######################
    
    if not only_single:
        if combined_fit:
            perfectClv = np.float(len(ydata[0][0]))
            perfectOn = np.float(len(ydata[0][1]))
            singleClv = 0.0
            singleOn = 0.0
            doubleClv = 0.0
            doubleOn = 0.0
            for i in range(len(xdata)):
                if len(xdata[i])==1:
                    if xdata[i][0] == 6 or xdata[i][0] == 7: #Here I eliminate some of the Boyle data points, comment out if not necessary
                        ydata[i][1] = []; yerr[i][1] = [];
                    singleClv += len(ydata[i][0])
                    singleOn += len(ydata[i][1])
                if len(xdata[i])==2:
                    if xdata[i][0] == 6 or xdata[i][0] == 7 or xdata[i][1] == 6 or xdata[i][1] == 7: #And here
                        ydata[i][1] = []; yerr[i][1] = [];
                    doubleClv += len(ydata[i][0])
                    doubleOn += len(ydata[i][1])

            chi_weights = [1/perfectClv,1/singleClv,1/doubleClv,1/perfectOn,1/singleOn,1/doubleOn]

        else:
            perfectClv = np.float(len(ydata[0]))
            singleClv = 0.0
            doubleClv = 0.0
            for i in range(len(xdata)):
                if len(xdata[i])==1:
                    singleClv += len(ydata[i])
                if len(xdata[i])==2:
                    doubleClv += len(ydata[i])

            chi_weights = [1/perfectClv,1/singleClv,1/doubleClv]
            
    else:
        if combined_fit:
            perfectClv = np.float(len(ydata[0][0]))
            perfectOn = np.float(len(ydata[0][1]))
            singleClv = 0.0
            singleOn = 0.0
            for i in range(len(xdata)):
                if len(xdata[i])==1:
                    if xdata[i][0]==1 or xdata[i][0]==3 or xdata[i][0]==4 or xdata[i][0]==6 or xdata[i][0]==7: #And here
                        ydata[i][1] = []
                        yerr[i][1] = []
                    singleClv += len(ydata[i][0])
                    singleOn += len(ydata[i][1])

            chi_weights = [1/perfectClv,1/singleClv,1/perfectOn,1/singleOn]
            
        else:
            perfectClv = np.float(len(ydata[0]))
            singleClv = 20.
            chi_weights = [1/perfectClv,1/singleClv,0.]
    
    
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
                                Tstart=80000.,             # infered from run on my computer/other runs on cluster
                                use_relative_steps=False,
                                delta=0.1,
                                tol=1E-5,
                                Tfinal=0.0,
                                potential_threshold = 300., #used when reanneal is True, above this threshold no solution is accepted
                                adjust_factor=1.1,
                                cooling_rate_high=0.95, #cooling rate for high T
                                cooling_rate_low=0.99, #cooling rate for low T
                                N_int=1000,
                                Ttransition=1000., #at this temp, change cooling rate
                                AR_low=40,
                                AR_high=60,
                                use_multiprocessing=True,
                                nprocs=nmbr_cores,
                                output_file_results = fit_result_file,
                                output_file_monitor = monitor_file,
                                output_file_init_monitor=init_monitor_file,
                                chi_weights=chi_weights, #automatically determined above
                                NMAC=False, #non-monotonic adaptive cooling, tested this to avoid local minima, but does not do much
                                reanneal=False, #reheating when in local minimum, set to False to do no reheating
                                combined_fit=combined_fit, #set in job file
                                random_step=False #every cycle takes steps not for all parameters, but only for randomly chosen set 
                                )

    t2 = time()

    print "Time elapsed: " + str((t2-t1)/(3600.)) + ' hrs'

    return






if __name__ == "__main__":
    main(sys.argv)


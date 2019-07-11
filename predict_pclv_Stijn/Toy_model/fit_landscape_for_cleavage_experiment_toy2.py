import numpy as np
import functools
import sys
import multiprocessing as mp

PATH_HPC05 = '/home/svandersmagt/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)
sys.path.append('../')
sys.path.append('../../code_Boyle')
import Nucleaseq_data_processing as processing
import calculate_cleavage_rate_toy as CRISPR
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
    use_cluster = bool(int(argv[7]))
    if use_cluster:
        path_to_dataClv= PATH_HPC05+'data_nucleaseq_Finkelsteinlab/targetE/'
        path_to_dataOn= path_to_dataClv
        nmbr_cores = 19
    else:
        # On my laptop use this line in stead:
        path_to_dataClv = '../../' + '/data_nucleaseq_Finkelsteinlab/targetE/'
        path_to_dataOn = '../../Data_Boyle/'
        nmbr_cores = 2

    # Simmulated Annealing
    combined_fit = bool(int(argv[8]))
    
    filename = argv[1]
    if combined_fit:
        model_ID =  argv[2] + '+' + argv[3]
    else:
        model_ID = argv[2]
    monitor_file = argv[4]
    fit_result_file = argv[5]
    init_monitor_file = argv[6]
    wa = True
    only_single = False
    log_boyle = True
    fill_in_boyle = True
    
    
    #gRNA_length = 20
    #fit_to_median = False   
    
    upper_bnd =      [5.0]  + [15.0]  + [7.0]  + [-1.0] + [3.0]*2  + [4.0]
    lower_bnd =      [-5.0] + [0.0]   + [3.0]  + [-3.0] + [-3.0]*2 + [-1.0]
    initial_guess =  [0.0]  + [0.0]   + [5.0]  + [-2.0] + [0.0]*2  + [2.0]
    

    ###########################
    # /* Objective function *\#
    ###########################
    KineticModel = functools.partial(CRISPR.calc_chi_squared,model_id=model_ID,log_on=log_boyle)


    #############################################
    # /* Preprocess the data from Boyle et al. *\#
    ##############################################
    if combined_fit:
        xdata, ydata, yerr = processing.prepare_multiprocessing_combined('1',filename,path_to_dataOn,path_to_dataClv,wa,only_single,log_boyle,fill_in_boyle)
    
    if not combined_fit:
        xdata, ydata, yerr = processing.prepare_multiprocessing_nucleaseq_log(filename,path_to_dataClv,True)
    #xdata, ydata, yerr = cr.create_fake_data()
    #print xdata, ydata, yerr
    # print ydata
    # print "test ... " + ' \n'
    # KineticModel(np.array(initial_guess),xdata,ydata,np.array([]),1.0)
    
    
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
                #if len(xdata[i])==1:
                #    if xdata[i][0]==6 or xdata[i][0]==7:
                #        ydata[i][1] = []
                #        yerr[i][1] = []
                if len(xdata[i]) > 0 and (xdata[i][0] == 2 or xdata[i][0] == 6 or xdata[i][0] == 7):
                    ydata[i][1] = []
                    yerr[i][1] = []
                if len(xdata[i])==1:
                    singleClv += len(ydata[i][0])
                    singleOn += len(ydata[i][1])
                if len(xdata[i])==2:
                    doubleClv += len(ydata[i][0])
                    doubleOn += len(ydata[i][1])

            chi_weights = [1/perfectClv,1/singleClv,1/doubleClv,1/perfectOn,1/singleOn,1/doubleOn]
            #chi_weights = [1.0,1.0,1.0,1.0,1.0,1.0] 
            #perfClv, sinClv, doubClv, perfOn, sinOn, doubOn

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
                if len(xdata[i]) > 0 and (xdata[i][0] == 6 or xdata[i][0] == 7):
                    ydata[i][1] = []
                    yerr[i][1] = []
                if len(xdata[i])==1:
                    #if xdata[i][0]==6 or xdata[i][0]==7:
                    #    ydata[i][1] = []
                    #    yerr[i][1] = []
                    singleClv += len(ydata[i][0])
                    singleOn += len(ydata[i][1])

            chi_weights = [1/perfectClv,1/singleClv,1/perfectOn,1/singleOn]
    
    
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
                                Tstart=20000.,             # infered from run on my computer/other runs on cluster
                                use_relative_steps=False,
                                delta=0.1,
                                tol=1E-5,
                                Tfinal=0.0,
                                potential_threshold = 1500.,
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
                                chi_weights=chi_weights,
                                NMAC=False, #non-monotonic adaptive cooling
                                reanneal=False, #reheating when in local minimum, set to False to do no reheating
                                combined_fit=combined_fit,
                                random_step=False
                                )

    t2 = time()

    print "Time elapsed: " + str((t2-t1)/(3600.)) + ' hrs'

    return






if __name__ == "__main__":
    main(sys.argv)

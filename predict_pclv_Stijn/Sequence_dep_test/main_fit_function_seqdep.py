import numpy as np
import functools
import sys
import multiprocessing as mp

PATH_HPC05 = '/home/svandersmagt/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)
sys.path.append(PATH_HPC05 + 'predict_pclv_Stijn')
sys.path.append('../../predict_pclv_Stijn')
sys.path.append('../../code_general')
import calculate_cleavage_rate as clv
import preprocessing as processing
import cleavage_rate as CRISPR
import SimulatedAnnealing_Nucleaseq_seqdep as SA

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
        #path_to_dataOn= path_to_dataClv
        nmbr_cores = 19
    else:
        # On my laptop use this line in stead:
        path_to_dataClv = '../../' + '/data_nucleaseq_Finkelsteinlab/'
        #path_to_dataOn = '../../Data_Boyle/'
        nmbr_cores = 2

    # Simmulated Annealing
    filename = argv[1]
    model_ID =  argv[2]
    monitor_file = argv[4]
    fit_result_file = argv[5]
    init_monitor_file = argv[6]
    
    gRNA_length = 20
    #fit_to_median = False   
    
    #indep
    upper_bnd = [10.]*40 + [6.] +[6.]
    lower_bnd = [-10.]*20 + [0.]*20 + [1.] + [3.]
    initial_guess = [0.]*20 + [5.]*20 + [3.] + [4.5]
    
    #v3
    #upper_bnd = [10.]*12 + [6.]*2
    #lower_bnd = [0.]*12 + [1.] + [3.]
    #initial_guess = [5.]*12 + [3.] + [4.5]
    
    #v2
    #upper_bnd = [10.0]*20 + [10.0]*16 + [4.0] + [6.0]
    #lower_bnd = [-10.0]*20 + [0.0]*16 + [1.0] + [3.0]
    #initial_guess =  [0.0]*20 + [5.0]*16 + [2.5] + [4.5]
    #initial_guess[20+1] = 1.5
    #initial_guess[20+4] = 1.5
    #initial_guess[20+11] = 1.5
    #initial_guess[20+14] = 1.5
    #upper_bnd[20+1] = 3.
    #upper_bnd[20+4] = 3.
    #upper_bnd[20+11] = 3.
    #upper_bnd[20+14] = 3.
    

    ###########################
    # /* Objective function *\#
    ###########################
    ##CHANGE THIS BACK
    KineticModel = functools.partial(clv.calc_chi_squared, chi_weights=[np.nan],
                        guide_length=gRNA_length,
                        model_id=model_ID)


    #############################################
    # /* Preprocess the data from Boyle et al. *\#
    ##############################################
    xdata1, ydata1, yerr1 = processing.prepare_multiprocessing_seq_indep(filename,path_to_dataClv)
    #xdata, ydata, yerr = cr.create_fake_data()
    #print xdata, ydata, yerr
    # print ydata
    # print "test ... " + ' \n'
    # KineticModel(np.array(initial_guess),xdata,ydata,np.array([]),1.0)
    xdata = []
    ydata = []
    yerr = []
    
    ##CHANGE THIS BACK
    for i in range(len(xdata1)):
        if len(xdata1[i])<2:
            xdata.append(xdata1[i])
            ydata.append(ydata1[i])
            yerr.append(yerr1[i])
    
    #chi_weights = [1.0,1.0,1.0,1.0,1.0,1.0] 
    #perfClv, sinClv, doubClv, perfOn, sinOn, doubOn 
    
    
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
                                Tstart=1000.,             # infered from run on my computer/other runs on cluster
                                use_relative_steps=False,
                                delta=1.0,
                                tol=1E-5,
                                Tfinal=0.,
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


import numpy as np
import functools
import sys
import multiprocessing as mp

PATH_HPC05 = '/home/svandersmagt/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)
import Nucleaseq_data_processing as processing
import calculate_cleavage_rate as CRISPR
import SimulatedAnnealing_Nucleaseq_parallel as SA
import create_fake_data as cr

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
        path_to_data= PATH_HPC05+'data_nucleaseq_Finkelsteinlab/targetE/'
        nmbr_cores = 19
    else:
        # On my laptop use this line in stead:
        path_to_data = '../' + '/data_nucleaseq_Finkelsteinlab/targetE/'
        nmbr_cores = 2

    # Simmulated Annealing
    filename = argv[1]
    model_ID =  argv[2]
    monitor_file = argv[3]
    fit_result_file = argv[4]
    init_monitor_file = argv[5]
    
    gRNA_length = 20
        
    
    upper_bnd = [15.0]*40 +  [7.0] *3
    lower_bnd = [-10.0]*20 + [0.0]*20 + [-5.0] *3
    initial_guess =  [0.0]*20 + [5.0]*20 + [0.0] *3
    

#    initial_guess = np.array([-0.91105443,  1.65099349, -6.28982814,  6.34251357,  0.24131213,  5.26027237,
#   1.05829753, -7.28159225, -2.86582133,  7.22669828,  6.12385666, -3.10445699,
#   8.8801378,  -3.97414648,  6.06113597, -8.51377994,  9.19863279,  6.59524463,
#   0.72173675, -8.47143605,  1.97290734,  6.28375984,  5.98784653,  4.02852934,
#   6.1261268,   0.55975488,  3.26946316,  4.09596255,  5.73483408,  0.37742956,
#   3.62945,     3.87178546,  3.87498714,  2.07444235,  1.45314807,  4.91403281,
#   2.8929454,   7.93602245,  1.64840129,  9.36968467, -2.03892706,  1.81355167,
#  -7.9])

    ###########################
    # /* Objective function *\#
    ###########################
    KineticModel = functools.partial(CRISPR.calc_chi_squared,
                        guide_length=gRNA_length,
                        model_id=model_ID)


    #############################################
    # /* Preprocess the data from Boyle et al. *\#
    ##############################################
    xdata, ydata, yerr = processing.prepare_multiprocessing_nucleaseq(filename,path_to_data)
    #xdata, ydata, yerr = cr.create_fake_data()
    #print xdata, ydata, yerr
   
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
                                Tstart=1000.,             # infered from run on my computer/other runs on cluster
                                use_relative_steps=False,
                                delta=0.1,
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


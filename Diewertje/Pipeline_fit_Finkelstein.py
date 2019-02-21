import numpy as np
import functools
import sys

PATH_HPC05 = '/home/mklein1/Energy_Landscape_dCas9/' #### Adjust this! Also in other files!
sys.path.append(PATH_HPC05)
import Prepare_data
import Chisq_Finkelstein as Chi
import Calculate_ABA_Finkelsteinlap as ABA
import SimulatedAnnealing_Finkelstein_parallel as SA


from time import time
'''
********************
Main script to be run on cluster for Simmulated Annealing fitting of energy landscape to data from Finkelstein.

Diewertje Dekker
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
        path_to_Boyle_data= PATH_HPC05+'/Data_Boyle/'
        nmbr_cores = 19
    else:
        # On my laptop use this line in stead:
        path_to_Boyle_data = '../' + '/Data_Boyle/'
        nmbr_cores = 1

    # Simmulated Annealing
    replica_ID = argv[1]
    model_ID =  argv[2]
    monitor_file = argv[3]
    fit_result_file = argv[4]
    init_monitor_file = argv[5]
    gRNA_length = 20


    upper_bnd =  [10.0] + [10.0]*40 +  [3.0] *2
    lower_bnd = [0.0] + [-10.0]*40 + [-7.0] *2
    initial_guess =  [5.0] + [0.0]*40 + [0.0] *2

    ###########################
    # /* Objective function *\#
    ###########################
    KineticModel = functools.partial(Chi.calc_Chi_square,
                        #guide_length=gRNA_length,
                        model_id=model_ID)
    # THIS ONE I NEEDED TO MAKE!
    
    concentrations=[0, 0.1, 0.3, 1, 3, 10, 30, 100, 300] # in nanoMolair
    reference=1 # in nanomolair
    #ontarget_ABA = functools.partial(ABA.calc_ABA(parameters, concentrations, reference, mismatch_positions=None, model_id = modelID, guide_length = gRNA_length, T=10*60), model_id=model_ID)

    #############################################
    # /* Preprocess the data from Finkelstein *\#
    #############################################
    path='../Data_ABA_Finkelsteinlab/'
    filename='cas9-target-e-replicate-1-delta-abas_Canonical_OT-r_0-20.csv'
    xdata, ydata, yerr = Prepare_data.Prepare_Cdata(path,filename)
    # xdata=MMpos, ydata=Delta ABA, yerr=Uncertainty


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
                                on_target_function=ABA.calc_ABA,
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

    print("Time elapsed: " + str((t2-t1)/(3600.)) + ' hrs')

    return


if __name__ == "__main__":
    main(sys.argv)

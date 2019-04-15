import numpy as np
import functools
import sys

PATH_HPC05 = '/home/dddekker/BEP' #### Adjust this! Also in other files!
#PATH_HPC05 = '/home/mklein1/Diewertje'
sys.path.append(PATH_HPC05)
import Prepare_data_Simple as Pre # for this rawABA data! for delta ABA use the normal one
import Chisq_Finkelstein as Chi
import SimulatedAnnealing_Finkelstein_parallel as SA
import Calculate_ABA_Finkelsteinlab_Diewertje as ABA

from time import time
'''
********************
Main script to be run on cluster for Simmulated Annealing fitting of energy landscape to data from Finkelstein.

Diewertje Dekker
********************

1. Load the datafiles from Finkelstein and preprocess it such that it allows for multiprocessing 
2. feed data and model into SA code 
'''


def main(argv):
    #################
    # /* Settings *\#
    #################
    # Loading the data
    use_cluster = bool(int(argv[5]))
    if use_cluster:
        path_to_data= PATH_HPC05 + '/Data_ABA_Finkelsteinlab/'
        nmbr_cores = 19
    else:
        # On my laptop use this line in stead:
        path_to_data = '../Data_ABA_Finkelsteinlab/' 
        nmbr_cores = 1

    # Simmulated Annealing
    model_ID =  argv[1] #'init_limit_general_energies_v2'
    monitor_file = argv[2] #'monitor.txt'
    fit_result_file = argv[3] #'fit_results.txt'
    init_monitor_file = argv[4] #'init_monitor.txt'
    gRNA_length = 20

    upper_bnd =  [10.0] + [10.0]*40 +  [3.0] *3 #*3
    lower_bnd = [0.0] + [-10.0]*40 + [-7.0] *3 #*3
    initial_guess = [5.0] + [3.0]*40 + [1.5] *3 #np.loadtxt('parameters.txt')

    ###########################
    # /* Objective function *\#
    ###########################
    concentrations = np.array([10,100])#[0.1, 0.3, 1, 3, 10, 30, 100, 300]) # in nanoMolair
    reference=1 # in nanomolair
    
    KineticModel = functools.partial(Chi.calc_Chi_square,model_id=model_ID, guide_length=gRNA_length,
                                     concentrations=concentrations,reference=reference)
     
    ###########################
    # /* Ontarget function *\#
    ###########################
    ONtarget=functools.partial(ABA.calc_ABA,concentrations=concentrations,reference=reference, 
                               mismatch_positions=[],model_id=model_ID, guide_length = gRNA_length, T=10*60)

    #############################################
    # /* Preprocess the data from Finkelstein *\#
    #############################################
    filename= 'TargetE-dCas9_AbsoluteABA_Canonical_OT-r_0-2.csv'
    xdata,ydata,yerr=Pre.Prepare_Cdata(path_to_data,filename) 
    # xdata=MMpos, ydata=ABA, yerr=Uncertainty
    
    #filename= 'cas9-target-e-replicate-1-delta-abas_Canonical_OT-r_0-2.csv'
    #xdata,ydata,yerr=Pre.Prepare_Cdata(path_to_data,filename)
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
                                    on_target_function=ONtarget,
                                    Tstart=100.,             # infered from run on my computer/other runs on cluster
                                    use_relative_steps=False,
                                    delta=1.0,
                                    tol=1E-5,
                                    Tfinal=0, #50
                                    adjust_factor=1.1,
                                    cooling_rate=0.99, #0.6
                                    N_int=1000, #10
                                    AR_low=40,
                                    AR_high=60,
                                    use_multiprocessing=True,
                                    nprocs=nmbr_cores,
                                    output_file_results = fit_result_file,
                                    output_file_monitor = monitor_file,
                                    output_file_init_monitor= init_monitor_file
                                  )

    t2 = time()

    print("Time elapsed: " + str((t2-t1)/(3600.)) + ' hrs')

    return


if __name__ == "__main__":
    main(sys.argv)


  ##############################################
  # /*   To run the code from the prompt     *\#
  ##############################################
# C:\Users\Diewertje\Documents\Year 3\BEP\Energy_Landscape_dCas9\Diewertje>python Pipeline_fit_Finkelstein.py 'init_limit_general_energies_v2' 'monitor.txt' 'fit_results.txt' 'init_monitor.txt' 0
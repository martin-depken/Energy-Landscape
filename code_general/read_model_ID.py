import numpy as np

choose_model = ['general_energies_rates',
                'general_energies',
                'init_limit_general_energies',
                'init_limit_fast_internal_general_energies',
                'constant_eps_I',
                'init_limit_lock_const_EpsI',
                'init_limit_two_drops_fixed_BP',
                'init_limit_two_drops',
                'init_limit_5EpsC_2EpsI']



def unpack_parameters(parameters, model_id='general_energies',guide_length=20):
    '''
    Use model ID to construct vector of epsilon values and forward rates.

    For every paramaterisation add a new case/ model_id
    :param parameters:
    :param model_id:
    :param guide_length:
    :return:
    '''

    epsilon = np.zeros(2 * guide_length + 1)
    forward_rates = np.ones(guide_length + 2)
    
    if model_id == 'Clv_Saturated_fixed_kf_general_energies_v2':
        if len(parameters)!=41:
            print 'Wrong number of parameters'
            return
        
        epsilon[0] = -100.0 #predefined epsilon PAM at saturation
        epsilon[1:] = parameters[:-1]
        
        rate_sol_to_PAM = 1000.0 #predefined at saturation
        rate_PAM_to_R1 = 200.0
        rate_internal = 200.0
        rate_clv = 10**parameters[-1]
        
        forward_rates = forward_rates * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = rate_clv
    
    elif model_id == 'Clv_Saturated_general_energies_v2':
        if len(parameters)!=42:
            print 'Wrong number of parameters'
            return
        
        epsilon[0] = -100.0 #predefined epsilon PAM at saturation
        epsilon[1:] = parameters[:-2]
        
        rate_sol_to_PAM = 1000.0 #predefined at saturation
        rate_PAM_to_R1 = 10**parameters[-2]
        rate_internal = 10**parameters[-2] #rate from PAM to R is equal to internal rate
        rate_clv = 10**parameters[-1]
        
        forward_rates = forward_rates * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = rate_clv

    elif model_id == 'Clv_init_limit_Saturated_general_energies_v2':
        if len(parameters)!=43:
            print 'Wrong number of parameters'
            return
        
        epsilon[0] = -100.0 #predefined epsilon PAM at saturation
        epsilon[1:] = parameters[:-3]
        
        rate_sol_to_PAM = 1000.0 #predefined at saturation
        rate_PAM_to_R1 = 10**parameters[-3]
        rate_internal = 10**parameters[-2]
        rate_clv = 10**parameters[-1]
        
        forward_rates = forward_rates * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = rate_clv
                

    elif model_id == 'Clv_init_limit_general_energies_v2':
        # General position dependency
        epsilon = parameters[:-4]

        rate_sol_to_PAM = 10**parameters[-4]
        rate_PAM_to_R1  = 10**parameters[-3]
        rate_internal = 10**parameters[-2]
        rate_clv = 10**parameters[-1]

        forward_rates = np.ones(guide_length + 2) * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = rate_clv

    
    elif model_id == 'general_energies_no_kPR':
        # ---- have the rate from PAM into R-loop the same as the forward rate within R-loop
        if len(parameters)!=43:
            print 'Wrong number of parameters'
            return
        # General position dependency
        epsilon = parameters[:-2]

        # --- rates: sol->PAM (concentration dependent), 1 constant forward rate for all remaining transitions
        rate_sol_to_PAM = 10**parameters[-2]
        rate_internal = 10**parameters[-1]

        forward_rates = np.ones(guide_length + 2) * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave


    elif model_id == 'landscape_lowest_chi_squared_fit_rates':
        # ---- fix the energies---
        # (copied from parameter file: '../data/25_10_2018/fit_25_10_2018_sim_22.txt') ----
        epsilon = np.array([  1.43364597e+00,  -2.51895658e+00,  -8.38107740e-01,
        -1.00837871e+00,  -3.89888343e+00,   4.98565931e+00,
        -2.24062010e+00,   1.75709991e+00,   1.48346110e+00,
        -2.56251518e+00,   4.76022290e+00,   1.66832631e+00,
        -4.41487326e-04,  -3.01917678e+00,   1.70186470e+00,
        -2.69692160e+00,   4.63508021e+00,   3.43845249e+00,
        -3.53360655e+00,   3.90785543e+00,   3.95624011e+00,
         8.41041112e+00,   3.52511767e+00,   6.47092824e+00,
         6.29617812e+00,   5.87466899e+00,   4.02069468e+00,
         6.97289538e+00,   5.39037459e+00,   6.53991724e+00,
         6.04624779e+00,   6.11140010e+00,   4.95893203e+00,
         5.40442705e+00,   5.69985755e+00,   5.12293027e+00,
         5.62074797e+00,   4.81777124e+00,   7.94515945e+00,
         9.77311952e+00,   6.84175107e+00])

       #  epsilon = np.array([ 1.31882561, -6.5880103 ,  1.59239502,  0.46021068, -2.49593644,
       #  0.09580053,  4.54430596, -3.37045113,  0.37192334,  1.02581499,
       #  4.12556609,  1.64960851, -0.03692466, -4.49653651,  4.39600456,
       # -3.57616013,  3.90152848,  3.48127153, -4.66585257,  1.77729046,
       #  8.90727104,  7.95522837,  4.24585854,  8.89394253,  8.89430408,
       #  5.15323997,  4.0149383 ,  6.61232836,  5.12389258,  7.22642299,
       #  6.06820965,  5.94807726,  4.90830081,  4.99741095,  6.38253949,
       #  5.87159526,  6.62698767,  5.87749165,  5.58373498,  9.01010833,
       #  4.79058499])

        # --- fit the timescales ----
        rate_sol_to_PAM = 10**parameters[0]
        rate_PAM_to_R1  = 10**parameters[1]
        rate_internal = 10**parameters[2]

        forward_rates = np.ones(guide_length + 2) * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = 0.0  # dCas9 does not cleave

    elif model_id == 'Boyle_median_landscape_fit_rates':

        # ---- fix the energies
        # (copied from parameter file: '../data/25_10_2018/median_landscape_Boyle_2Dgrid.txt' on 04/12/2018) ----
        epsilon = np.array([ 1.37168412, -4.15848621, -1.9680678 , -0.88508854, -0.70265355,
        2.00690677,  0.44846725, -0.98458337,  0.69054452,  1.39284825,
        3.99499433,  1.64210983,  0.01835355, -3.80696161,  2.1347165 ,
       -0.80030528,  1.97963966,  3.08286715, -0.99001786,  2.70016623,
        3.62856361,  7.15113618,  3.43168686,  8.05355657,  7.47981321,
        5.62454054,  3.95761041,  6.69713956,  5.23678395,  7.46232812,
        6.10345114,  6.1114001 ,  4.97393321,  5.29781208,  6.10472344,
        5.43762448,  5.02455717,  4.3662047 ,  3.43647645,  7.07048864,
        5.22717571])

        # --- fit the timescales ----
        rate_sol_to_PAM = 10**parameters[0]
        rate_PAM_to_R1  = 10**parameters[1]
        rate_internal = 10**parameters[2]

        forward_rates = np.ones(guide_length + 2) * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = 0.0  # dCas9 does not cleave


    elif model_id == 'general_energies_rates':
        epsilon = parameters[:(2*guide_length+1)]
        forward_rates = np.ones(guide_length + 2)
        forward_rates[:-1] = 10**np.array(parameters[(2*guide_length+1):])
        forward_rates[-1] = 0.0 # dCas9 does not cleave

    elif model_id == 'general_energies':
        # General position dependency + minimal amount of rates
        epsilon = parameters[:-2]
        forward_rates = np.ones(guide_length + 2) * parameters[-2] #internal rates
        forward_rates[0] = parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave

    elif model_id == 'general_energies_v2':
        # General position dependency + minimal amount of rates
        epsilon = parameters[:-2]
        forward_rates = np.ones(guide_length + 2) * 10**parameters[-2] #internal rates
        forward_rates[0] = 10**parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave

    elif model_id == 'init_limit_general_energies_v0':
        # General position dependency
        epsilon = parameters[:-3]
        forward_rates = np.ones(guide_length + 2) * parameters[-2] #internal rates
        forward_rates[0] = parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave
        # first rate from PAM into R-loop is less or equal to other internal rates
        forward_rates[1] = np.exp(-parameters[-3])* parameters[-2]


    elif model_id == 'init_limit_general_energies':
        # General position dependency
        epsilon = parameters[:-3]
        forward_rates = np.ones(guide_length + 2) * 10**parameters[-2] #internal rates
        forward_rates[0] = 10**parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave
        # first rate from PAM into R-loop is less or equal to other internal rates
        forward_rates[1] = np.exp(-parameters[-3])* 10**parameters[-2]

    elif model_id == 'init_limit_general_energies_v2':
        # General position dependency
        epsilon = parameters[:-3]

        rate_sol_to_PAM = 10**parameters[-3]
        rate_PAM_to_R1  = 10**parameters[-2]
        rate_internal = 10**parameters[-1]

        forward_rates = np.ones(guide_length + 2) * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = 0.0  # dCas9 does not cleave


    elif model_id == 'init_limit_fast_internal_general_energies':
        # General position dependency for energies
        epsilon = parameters[:-3]

        internal_rates = 10**parameters[-1]  #  directly set the rate (the order of magnitude)
        rate_sol_to_PAM = 10**parameters[-3]  # directly set a rate (the order of magnitude)
        rate_PAM_to_R1 = np.exp(-parameters[-2]+ epsilon[1]) * internal_rates   # through placement of the transition state above R1's energy.

        forward_rates = np.ones(guide_length + 2)* internal_rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = 0.0





    elif model_id == 'constant_eps_I':
        # General position dependency for matches, constant mismatch penalty
        epsPAM = parameters[0]
        epsilonC = parameters[1:(guide_length+1)]
        epsilonI = parameters[guide_length+1]

        epsilon[0] = epsPAM
        epsilon[1:(guide_length+1)] = epsilonC
        epsilon[(guide_length+1):] = epsilonI

        forward_rates = np.ones(guide_length + 2) * parameters[-2] #internal rates
        forward_rates[0] = parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave

    elif model_id == 'init_limit_lock_const_EpsI':
        e_PAM = parameters[0]
        ec_1 = parameters[1]
        ec_first = parameters[2]
        ec_second = parameters[3]
        e_I = parameters[4]
        x = parameters[5]
        k_PAM = parameters[6]
        E_barr = parameters[7]
        k = parameters[8]
        k_1 = k*np.exp(-E_barr)


        epsilon[0] = e_PAM
        epsilon[1] = ec_1
        epsilon[2:x + 1] = ec_first
        epsilon[x + 1:guide_length + 1] = ec_second
        epsilon[guide_length + 1:] = e_I

        forward_rates = np.ones(guide_length + 2) * k  # internal rates
        forward_rates[0] = k_PAM
        forward_rates[1] = k_1
        forward_rates[-1] = 0.0

    elif model_id == 'init_limit_two_drops_fixed_BP':
        pos1 = 10
        pos2 = 18

        e_PAM = parameters[0]
        ec_1 = parameters[1]
        ec_first = parameters[2]
        drop1 = parameters[3]
        ec_second = parameters[4]
        drop2 = parameters[5]
        ec_third = parameters[6]
        e_I = parameters[7]
        k_PAM = parameters[8]
        E_barr = parameters[9]
        k = parameters[10]
        k_1 = k*np.exp(-E_barr)

        epsilon[0] = e_PAM
        epsilon[1] = ec_1
        epsilon[2:pos1 + 1] = ec_first
        epsilon[pos1+1] = drop1
        epsilon[pos1 + 2:pos2 + 1] = ec_second
        epsilon[pos2+1]= drop2
        epsilon[pos2+2:guide_length+1] = ec_third
        epsilon[guide_length+1:] = e_I

        forward_rates = np.ones(guide_length + 2) * k  # internal rates
        forward_rates[0] = k_PAM
        forward_rates[1] = k_1
        forward_rates[-1] = 0.0

    elif model_id == 'init_limit_two_drops':

        e_PAM = parameters[0]
        ec_1 = parameters[1]
        ec_first = parameters[2]
        drop1 = parameters[3]
        ec_second = parameters[4]
        drop2 = parameters[5]
        ec_third = parameters[6]
        pos1 = parameters[7]
        pos2 = parameters[8]
        e_I = parameters[9]
        k_PAM = parameters[10]
        E_barr = parameters[11]
        k = parameters[12]
        k_1 = k * np.exp(-E_barr)

        epsilon[0] = e_PAM
        epsilon[1] = ec_1
        epsilon[2:pos1 + 1] = ec_first
        epsilon[pos1 + 1] = drop1
        epsilon[pos1 + 2:pos2 + 1] = ec_second
        epsilon[pos2 + 1] = drop2
        epsilon[pos2 + 2:guide_length + 1] = ec_third
        epsilon[guide_length + 1:] = e_I

        forward_rates = np.ones(guide_length + 2) * k  # internal rates
        forward_rates[0] = k_PAM
        forward_rates[1] = k_1
        forward_rates[-1] = 0.0

    elif model_id == 'init_limit_5EpsC_2EpsI':
        e_PAM = parameters[0]
        ec_1 = parameters[1]
        ec_2 = parameters[2]
        ec_3 = parameters[3]
        ec_4 = parameters[4]
        ec_5 = parameters[5]
        eI_1 = parameters[6]
        eI_2 = parameters[7]
        bp2 = parameters[8]
        bp3 = parameters[9]
        bp4 = parameters[10]
        bpI = parameters[11]
        k_PAM = parameters[12]
        E_barr = parameters[13]
        k = parameters[14]
        k_1 = k * np.exp(-E_barr)

        epsilon[0] = e_PAM
        epsilon[1] = ec_1
        epsilon[2:bp2 + 1] = ec_2
        epsilon[bp2+1:bp3 + 1] = ec_3
        epsilon[bp3 + 1:bp4 + 1] = ec_4
        epsilon[bp4 + 1:guide_length + 1] = ec_5
        epsilon[guide_length + 1:guide_length + bpI + 1] = eI_1
        epsilon[guide_length + bpI + 1:] = eI_2

        forward_rates = np.ones(guide_length + 2) * k  # internal rates
        forward_rates[0] = k_PAM
        forward_rates[1] = k_1
        forward_rates[-1] = 0.0
    else:
        print 'Watch out! Non-existing model-ID..'
        return



    return epsilon, forward_rates
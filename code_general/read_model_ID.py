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

    if model_id == 'general_energies_rates':
        epsilon = parameters[:(2*guide_length+1)]
        forward_rates = np.ones(guide_length + 2)
        forward_rates[:-1] = 10**np.array(parameters[(2*guide_length+1):])
        forward_rates[-1] = 0.0 # dCas9 does not cleave

    if model_id == 'general_energies':
        # General position dependency + minimal amount of rates
        epsilon = parameters[:-2]
        forward_rates = np.ones(guide_length + 2) * parameters[-2] #internal rates
        forward_rates[0] = parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave

    if model_id == 'general_energies_v2':
        # General position dependency + minimal amount of rates
        epsilon = parameters[:-2]
        forward_rates = np.ones(guide_length + 2) * 10**parameters[-2] #internal rates
        forward_rates[0] = 10**parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave

    if model_id == 'init_limit_general_energies_v0':
        # General position dependency
        epsilon = parameters[:-3]
        forward_rates = np.ones(guide_length + 2) * parameters[-2] #internal rates
        forward_rates[0] = parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave
        # first rate from PAM into R-loop is less or equal to other internal rates
        forward_rates[1] = np.exp(-parameters[-3])* parameters[-2]


    if model_id == 'init_limit_general_energies':
        # General position dependency
        epsilon = parameters[:-3]
        forward_rates = np.ones(guide_length + 2) * 10**parameters[-2] #internal rates
        forward_rates[0] = 10**parameters[-1]  # from solution to PAM
        forward_rates[-1] = 0.0  # dCas9 does not cleave
        # first rate from PAM into R-loop is less or equal to other internal rates
        forward_rates[1] = np.exp(-parameters[-3])* 10**parameters[-2]

    if model_id == 'init_limit_general_energies_v2':
        # General position dependency
        epsilon = parameters[:-3]

        rate_sol_to_PAM = 10**parameters[-3]
        rate_PAM_to_R1  = 10**parameters[-2]
        rate_internal = 10**parameters[-1]

        forward_rates = np.ones(guide_length + 2) * rate_internal #internal rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = 0.0  # dCas9 does not cleave


    if model_id == 'init_limit_fast_internal_general_energies':
        # General position dependency for energies
        epsilon = parameters[:-3]

        internal_rates = 10**parameters[-1]  #  directly set the rate (the order of magnitude)
        rate_sol_to_PAM = 10**parameters[-3]  # directly set a rate (the order of magnitude)
        rate_PAM_to_R1 = np.exp(-parameters[-2]+ epsilon[1]) * internal_rates   # through placement of the transition state above R1's energy.

        forward_rates = np.ones(guide_length + 2)* internal_rates
        forward_rates[0] = rate_sol_to_PAM
        forward_rates[1] = rate_PAM_to_R1
        forward_rates[-1] = 0.0





    if model_id == 'constant_eps_I':
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

    if model_id == 'init_limit_lock_const_EpsI':
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

    if model_id == 'init_limit_two_drops_fixed_BP':
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

    if model_id == 'init_limit_two_drops':

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

    if model_id == 'init_limit_5EpsC_2EpsI':
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

    return epsilon, forward_rates
import numpy as np



def unpack_parameters(parameters, model_id,guide_length=20):
    '''
    Use model ID to construct vector of epsilon values and forward rates.

    For every paramaterisation add a new case/ model_id
    :param parameters:
    :param model_id:
    :param guide_length:
    :return:
    '''
    
    #si = sequence independent part, fit 18-7 nr 17
    epsilon_si = np.array([-1.00000000e+02, -6.33371589e+00, -1.36490455e+00,  4.44466676e+00,
                     -1.53708382e+00,  8.52684738e-01, -5.59779326e-02, -2.48652563e+00,
                      1.43245925e+00,  3.75567194e+00,  1.15654484e+00,  1.02454421e+00,
                      6.19640712e-01, -3.13409723e+00,  1.64314626e+00,  1.07824129e+00,
                     -2.15390976e+00,  1.04251947e+00,  3.02170067e+00, -3.99089331e-01,
                      6.78994364e+00,  5.65333843e+00,  4.11041237e+00,  6.48240764e+00,
                      6.97672083e+00,  6.26648141e+00,  7.38856778e+00,  6.89948764e+00,
                      6.22048036e+00,  8.99420177e+00,  7.25361518e+00,  7.40355908e+00,
                      7.02417975e+00,  7.73956089e+00,  7.88443131e+00,  7.64480966e+00,
                      6.35536811e+00,  5.13351785e+00,  4.24795058e+00,  5.83044258e+00,
                      2.41878155e+00])
    rates_si = np.array([1000.,          641.2890381,   641.2890381,   641.2890381,   641.2890381,
                  641.2890381,   641.2890381,   641.2890381,   641.2890381,   641.2890381,
                  641.2890381,   641.2890381,   641.2890381,   641.2890381,   641.2890381,
                  641.2890381,   641.2890381,   641.2890381,   641.2890381,   641.2890381,
                  641.2890381,   2.39286365])
    
    epsilon_si_aba = np.array([2.73766502381, -6.33371589e+00, -1.36490455e+00,  4.44466676e+00,
                     -1.53708382e+00,  8.52684738e-01, -5.59779326e-02, -2.48652563e+00,
                      1.43245925e+00,  3.75567194e+00,  1.15654484e+00,  1.02454421e+00,
                      6.19640712e-01, -3.13409723e+00,  1.64314626e+00,  1.07824129e+00,
                     -2.15390976e+00,  1.04251947e+00,  3.02170067e+00, -3.99089331e-01,
                      6.78994364e+00,  5.65333843e+00,  4.11041237e+00,  6.48240764e+00,
                      6.97672083e+00,  6.26648141e+00,  7.38856778e+00,  6.89948764e+00,
                      6.22048036e+00,  8.99420177e+00,  7.25361518e+00,  7.40355908e+00,
                      7.02417975e+00,  7.73956089e+00,  7.88443131e+00,  7.64480966e+00,
                      6.35536811e+00,  5.13351785e+00,  4.24795058e+00,  5.83044258e+00,
                      2.41878155e+00])
    rates_si_aba = np.array([8.58596308e-04,    641.2890381,   641.2890381,   641.2890381,   641.2890381,
                  641.2890381,   641.2890381,   641.2890381,   641.2890381,   641.2890381,
                  641.2890381,   641.2890381,   641.2890381,   641.2890381,   641.2890381,
                  641.2890381,   641.2890381,   641.2890381,   641.2890381,   641.2890381,
                  641.2890381,   0.0])
    
    if model_id == 'seq_dep_finetuning_direct_pairs' or model_id == 'seq_dep_finetuning_direct_pairs_aba' or model_id == 'seq_dep_finetuning_direct_pairs_constant_ei' or model_id == 'seq_dep_finetuning_direct_pairs_constant_ei_aba': #finetuning of fit 18-7 nr 17
        if len(parameters) != 12:
            print 'Wrong number of parameters'
            return
        
        epsilon_i_sd = np.zeros(33) #sd = sequence dependent
        
        epsilon_i_sd[0]  = parameters[0]  #A-AA
        epsilon_i_sd[1]  = parameters[0]  #U-AA
        epsilon_i_sd[2]  = parameters[0]  #G-AA
        epsilon_i_sd[3]  = parameters[0]  #C-AA
        epsilon_i_sd[4]  = parameters[1]  #A-AG
        epsilon_i_sd[5]  = parameters[1]  #U-AG
        epsilon_i_sd[6]  = parameters[1]  #G-AG
        epsilon_i_sd[7]  = parameters[1]  #C-AG
        epsilon_i_sd[8]  = parameters[2]  #A-AC
        epsilon_i_sd[9]  = parameters[2]  #U-AC
        epsilon_i_sd[10] = parameters[2]  #G-AC
        epsilon_i_sd[11] = parameters[2]  #C-AC
        epsilon_i_sd[12] = parameters[3]  #A-UT
        epsilon_i_sd[13] = parameters[3]  #G-UT
        epsilon_i_sd[14] = parameters[4]  #A-UG
        epsilon_i_sd[15] = parameters[4]  #G-UG
        epsilon_i_sd[16] = parameters[5]  #A-UC
        epsilon_i_sd[17] = parameters[5]  #G-UC
        epsilon_i_sd[18] = parameters[6]  #A-GA
        epsilon_i_sd[19] = parameters[6]  #C-GA
        epsilon_i_sd[20] = parameters[7]  #A-GT
        epsilon_i_sd[21] = parameters[7]  #C-GT
        epsilon_i_sd[22] = parameters[8]  #A-GG
        epsilon_i_sd[23] = parameters[8]  #C-GG
        epsilon_i_sd[24] = parameters[9]  #A-CA
        epsilon_i_sd[25] = parameters[9]  #G-CA
        epsilon_i_sd[26] = parameters[9]  #PAM-CA
        epsilon_i_sd[27] = parameters[10] #A-CT
        epsilon_i_sd[28] = parameters[10] #G-CT
        epsilon_i_sd[29] = parameters[10] #PAM-CT
        epsilon_i_sd[30] = parameters[11] #A-CC
        epsilon_i_sd[31] = parameters[11] #G-CC
        epsilon_i_sd[32] = parameters[11] #PAM-CC
        
        if model_id == 'seq_dep_finetuning_direct_pairs':
            return epsilon_si,rates_si,epsilon_i_sd
        elif model_id == 'seq_dep_finetuning_direct_pairs_constant_ei':
            epsilon_si[-20:] = 0.
            return epsilon_si,rates_si,epsilon_i_sd
        elif model_id == 'seq_dep_finetuning_direct_pairs_constant_ei_aba':
            epsilon_si_aba[-20:] = 0.
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
        else:
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
    
    
    elif model_id == 'seq_dep_finetuning_nearest_neighbour' or model_id == 'seq_dep_finetuning_nearest_neighbour_aba' or model_id == 'seq_dep_finetuning_nearest_neighbour_constant_ei' or model_id == 'seq_dep_finetuning_nearest_neighbour_constant_ei_aba': #finetuning of fit 18-7 nr 17
        if len(parameters) != 33:
            print 'Wrong number of parameters'
            return
        
        epsilon_i_sd = np.zeros(33) #sd = sequence dependent
        
        epsilon_i_sd[0]  = parameters[0]  #A-AA
        epsilon_i_sd[1]  = parameters[1]  #U-AA
        epsilon_i_sd[2]  = parameters[2]  #G-AA
        epsilon_i_sd[3]  = parameters[3]  #C-AA
        epsilon_i_sd[4]  = parameters[4]  #A-AG
        epsilon_i_sd[5]  = parameters[5]  #U-AG
        epsilon_i_sd[6]  = parameters[6]  #G-AG
        epsilon_i_sd[7]  = parameters[7]  #C-AG
        epsilon_i_sd[8]  = parameters[8]  #A-AC
        epsilon_i_sd[9]  = parameters[9]  #U-AC
        epsilon_i_sd[10] = parameters[10] #G-AC
        epsilon_i_sd[11] = parameters[11] #C-AC
        epsilon_i_sd[12] = parameters[12] #A-UT
        epsilon_i_sd[13] = parameters[13] #G-UT
        epsilon_i_sd[14] = parameters[14] #A-UG
        epsilon_i_sd[15] = parameters[15] #G-UG
        epsilon_i_sd[16] = parameters[16] #A-UC
        epsilon_i_sd[17] = parameters[17] #G-UC
        epsilon_i_sd[18] = parameters[18] #A-GA
        epsilon_i_sd[19] = parameters[19] #C-GA
        epsilon_i_sd[20] = parameters[20] #A-GT
        epsilon_i_sd[21] = parameters[21] #C-GT
        epsilon_i_sd[22] = parameters[22] #A-GG
        epsilon_i_sd[23] = parameters[23] #C-GG
        epsilon_i_sd[24] = parameters[24] #A-CA
        epsilon_i_sd[25] = parameters[25] #G-CA
        epsilon_i_sd[26] = parameters[26] #PAM-CA
        epsilon_i_sd[27] = parameters[27] #A-CT
        epsilon_i_sd[28] = parameters[28] #G-CT
        epsilon_i_sd[29] = parameters[29] #PAM-CT
        epsilon_i_sd[30] = parameters[30] #A-CC
        epsilon_i_sd[31] = parameters[31] #G-CC
        epsilon_i_sd[32] = parameters[32] #PAM-CT
        
        if model_id == 'seq_dep_finetuning_nearest_neighbour':
            return epsilon_si,rates_si,epsilon_i_sd
        elif model_id == 'seq_dep_finetuning_nearest_neighbour_constant_ei':
            epsilon_si[-20:] = 0.
            return epsilon_si,rates_si,epsilon_i_sd
        elif model_id == 'seq_dep_finetuning_nearest_neighbour_constant_ei_aba':
            epsilon_si_aba[-20:] = 0.
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
        else:
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
    
    
    elif model_id == model_id == 'seq_dep_finetuning_direct_pairs_free_ei' or model_id == 'seq_dep_finetuning_direct_pairs_free_ei_aba':
        if len(parameters) != 32:
            print 'Wrong number of parameters'
            return
        
        epsilon_i_sd = np.zeros(33) #sd = sequence dependent
        
        epsilon_i_sd[0]  = parameters[0]  #A-AA
        epsilon_i_sd[1]  = parameters[0]  #U-AA
        epsilon_i_sd[2]  = parameters[0]  #G-AA
        epsilon_i_sd[3]  = parameters[0]  #C-AA
        epsilon_i_sd[4]  = parameters[1]  #A-AG
        epsilon_i_sd[5]  = parameters[1]  #U-AG
        epsilon_i_sd[6]  = parameters[1]  #G-AG
        epsilon_i_sd[7]  = parameters[1]  #C-AG
        epsilon_i_sd[8]  = parameters[2]  #A-AC
        epsilon_i_sd[9]  = parameters[2]  #U-AC
        epsilon_i_sd[10] = parameters[2]  #G-AC
        epsilon_i_sd[11] = parameters[2]  #C-AC
        epsilon_i_sd[12] = parameters[3]  #A-UT
        epsilon_i_sd[13] = parameters[3]  #G-UT
        epsilon_i_sd[14] = parameters[4]  #A-UG
        epsilon_i_sd[15] = parameters[4]  #G-UG
        epsilon_i_sd[16] = parameters[5]  #A-UC
        epsilon_i_sd[17] = parameters[5]  #G-UC
        epsilon_i_sd[18] = parameters[6]  #A-GA
        epsilon_i_sd[19] = parameters[6]  #C-GA
        epsilon_i_sd[20] = parameters[7]  #A-GT
        epsilon_i_sd[21] = parameters[7]  #C-GT
        epsilon_i_sd[22] = parameters[8]  #A-GG
        epsilon_i_sd[23] = parameters[8]  #C-GG
        epsilon_i_sd[24] = parameters[9]  #A-CA
        epsilon_i_sd[25] = parameters[9]  #G-CA
        epsilon_i_sd[26] = parameters[9]  #PAM-CA
        epsilon_i_sd[27] = parameters[10] #A-CT
        epsilon_i_sd[28] = parameters[10] #G-CT
        epsilon_i_sd[29] = parameters[10] #PAM-CT
        epsilon_i_sd[30] = parameters[11] #A-CC
        epsilon_i_sd[31] = parameters[11] #G-CC
        epsilon_i_sd[32] = parameters[11] #PAM-CC
        
        epsilon_si[-20:] = parameters[12:]
        epsilon_si_aba[-20:] = parameters[12:]
        
        if model_id == 'seq_dep_finetuning_direct_pairs_free_ei':
            return epsilon_si,rates_si,epsilon_i_sd
        else:
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
        
        
    elif model_id == 'seq_dep_finetuning_nearest_neighbour_free_ei' or model_id == 'seq_dep_finetuning_nearest_neighbour_free_ei_aba':
        if len(parameters) != 53:
            print 'Wrong number of parameters'
            return
        
        epsilon_i_sd = np.zeros(33) #sd = sequence dependent
        
        epsilon_i_sd[0]  = parameters[0]  #A-AA
        epsilon_i_sd[1]  = parameters[1]  #U-AA
        epsilon_i_sd[2]  = parameters[2]  #G-AA
        epsilon_i_sd[3]  = parameters[3]  #C-AA
        epsilon_i_sd[4]  = parameters[4]  #A-AG
        epsilon_i_sd[5]  = parameters[5]  #U-AG
        epsilon_i_sd[6]  = parameters[6]  #G-AG
        epsilon_i_sd[7]  = parameters[7]  #C-AG
        epsilon_i_sd[8]  = parameters[8]  #A-AC
        epsilon_i_sd[9]  = parameters[9]  #U-AC
        epsilon_i_sd[10] = parameters[10] #G-AC
        epsilon_i_sd[11] = parameters[11] #C-AC
        epsilon_i_sd[12] = parameters[12] #A-UT
        epsilon_i_sd[13] = parameters[13] #G-UT
        epsilon_i_sd[14] = parameters[14] #A-UG
        epsilon_i_sd[15] = parameters[15] #G-UG
        epsilon_i_sd[16] = parameters[16] #A-UC
        epsilon_i_sd[17] = parameters[17] #G-UC
        epsilon_i_sd[18] = parameters[18] #A-GA
        epsilon_i_sd[19] = parameters[19] #C-GA
        epsilon_i_sd[20] = parameters[20] #A-GT
        epsilon_i_sd[21] = parameters[21] #C-GT
        epsilon_i_sd[22] = parameters[22] #A-GG
        epsilon_i_sd[23] = parameters[23] #C-GG
        epsilon_i_sd[24] = parameters[24] #A-CA
        epsilon_i_sd[25] = parameters[25] #G-CA
        epsilon_i_sd[26] = parameters[26] #PAM-CA
        epsilon_i_sd[27] = parameters[27] #A-CT
        epsilon_i_sd[28] = parameters[28] #G-CT
        epsilon_i_sd[29] = parameters[29] #PAM-CT
        epsilon_i_sd[30] = parameters[30] #A-CC
        epsilon_i_sd[31] = parameters[31] #G-CC
        epsilon_i_sd[32] = parameters[32] #PAM-CT
        
        epsilon_si[-20:] = parameters[33:]
        epsilon_si_aba[-20:] = parameters[33:]
        
        if model_id == 'seq_dep_finetuning_nearest_neighbour_free_ei':
            return epsilon_si,rates_si,epsilon_i_sd
        else:
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
        
        
    elif model_id == 'seq_dep_finetuning_nearest_neighbour_after' or model_id == 'seq_dep_finetuning_nearest_neighbour_after_aba' or model_id == 'seq_dep_finetuning_nearest_neighbour_after_constant_ei' or model_id == 'seq_dep_finetuning_nearest_neighbour_after_constant_ei_aba': #finetuning of fit 18-7 nr 17
        if len(parameters) != 30:
            print 'Wrong number of parameters'
            return
        
        epsilon_i_sd = np.zeros(30) #sd = sequence dependent
        
        epsilon_i_sd[0]  = parameters[0]  #AA-A
        epsilon_i_sd[1]  = parameters[1]  #AA-U
        epsilon_i_sd[2]  = parameters[2]  #AA-G
        epsilon_i_sd[3]  = parameters[3]  #AA-C
        epsilon_i_sd[4]  = parameters[4]  #AG-A
        epsilon_i_sd[5]  = parameters[5]  #AG-U
        epsilon_i_sd[6]  = parameters[6]  #AG-G
        epsilon_i_sd[7]  = parameters[7]  #AG-C
        epsilon_i_sd[8]  = parameters[8]  #AC-A
        epsilon_i_sd[9]  = parameters[9]  #AC-U
        epsilon_i_sd[10] = parameters[10] #AC-G
        epsilon_i_sd[11] = parameters[11] #AC-C
        epsilon_i_sd[12] = parameters[12] #UT-A
        epsilon_i_sd[13] = parameters[13] #UG-A
        epsilon_i_sd[14] = parameters[14] #UC-A
        epsilon_i_sd[15] = parameters[15] #GA-A
        epsilon_i_sd[16] = parameters[16] #GA-U
        epsilon_i_sd[17] = parameters[17] #GA-C
        epsilon_i_sd[18] = parameters[18] #GT-A
        epsilon_i_sd[19] = parameters[19] #GT-U
        epsilon_i_sd[20] = parameters[20] #GT-C
        epsilon_i_sd[21] = parameters[21] #GG-A
        epsilon_i_sd[22] = parameters[22] #GG-U
        epsilon_i_sd[23] = parameters[23] #GG-C
        epsilon_i_sd[24] = parameters[24] #CA-A
        epsilon_i_sd[25] = parameters[25] #CA-G
        epsilon_i_sd[26] = parameters[26] #CT-A
        epsilon_i_sd[27] = parameters[27] #CT-G
        epsilon_i_sd[28] = parameters[28] #CC-A
        epsilon_i_sd[29] = parameters[29] #CC-G

        
        if model_id == 'seq_dep_finetuning_nearest_neighbour_after':
            return epsilon_si,rates_si,epsilon_i_sd
        elif model_id == 'seq_dep_finetuning_nearest_neighbour_after_constant_ei':
            epsilon_si[-20:] = 0.
            return epsilon_si,rates_si,epsilon_i_sd
        elif model_id == 'seq_dep_finetuning_nearest_neighbour_after_constant_ei_aba':
            epsilon_si_aba[-20:] = 0.
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
        else:
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
        
        
    elif model_id == 'seq_dep_finetuning_nearest_neighbour_both' or model_id == 'seq_dep_finetuning_nearest_neighbour_both_aba' or model_id == 'seq_dep_finetuning_nearest_neighbour_both_constant_ei' or model_id == 'seq_dep_finetuning_nearest_neighbour_both_constant_ei_aba': #finetuning of fit 18-7 nr 17
        if len(parameters) != 45:
            print 'Wrong number of parameters'
            return
        
        epsilon_i_sd = np.zeros(45) #sd = sequence dependent
        
        [
            'A-UC-A','G-UC-A','A-GA-A','A-GA-U','C-GA-C','A-GT-A','A-GT-U','C-GT-C','A-GG-A','A-GG-U','C-GG-C','A-CA-G',
            'G-CA-A','P-CA-G','A-CT-G','G-CT-A','P-CT-G','A-CC-G','G-CC-A','P-CC-G']

        epsilon_i_sd[0]  = parameters[27] #A-AA-A        
        epsilon_i_sd[1]  = parameters[30] #A-AA-U
        epsilon_i_sd[2]  = parameters[21] #U-AA-G
        epsilon_i_sd[3]  = parameters[36] #U-AA-C
        epsilon_i_sd[4]  = parameters[24] #G-AA-A
        epsilon_i_sd[5]  = parameters[12] #G-AA-G
        epsilon_i_sd[6]  = parameters[6]  #C-AA-G
        epsilon_i_sd[7]  = parameters[28] #A-AG-A
        epsilon_i_sd[8]  = parameters[31] #A-AG-U
        epsilon_i_sd[9]  = parameters[22] #U-AG-G
        epsilon_i_sd[10] = parameters[37] #U-AG-C
        epsilon_i_sd[11] = parameters[25] #G-AG-A
        epsilon_i_sd[12] = parameters[13] #G-AG-G
        epsilon_i_sd[13] = parameters[7]  #C-AG-G
        epsilon_i_sd[14] = parameters[29] #A-AC-A
        epsilon_i_sd[15] = parameters[32] #A-AC-U
        epsilon_i_sd[16] = parameters[23] #U-AC-G
        epsilon_i_sd[17] = parameters[38] #U-AC-C
        epsilon_i_sd[18] = parameters[26] #G-AC-A
        epsilon_i_sd[19] = parameters[14] #G-AC-G
        epsilon_i_sd[20] = parameters[8]  #C-AC-G
        epsilon_i_sd[21] = parameters[33] #A-UT-A
        epsilon_i_sd[22] = parameters[18] #G-UT-A
        epsilon_i_sd[23] = parameters[34] #A-UG-A
        epsilon_i_sd[24] = parameters[19] #G-UG-A
        epsilon_i_sd[25] = parameters[35] #A-UC-A
        epsilon_i_sd[26] = parameters[20] #G-UC-A
        epsilon_i_sd[27] = parameters[9]  #A-GA-A
        epsilon_i_sd[28] = parameters[15] #A-GA-U
        epsilon_i_sd[29] = parameters[0]  #C-GA-C
        epsilon_i_sd[30] = parameters[10] #A-GT-A
        epsilon_i_sd[31] = parameters[16] #A-GT-U
        epsilon_i_sd[32] = parameters[1]  #C-GT-C
        epsilon_i_sd[33] = parameters[11] #A-GG-A
        epsilon_i_sd[34] = parameters[17] #A-GG-U
        epsilon_i_sd[35] = parameters[2]  #C-GG-C
        epsilon_i_sd[36] = parameters[39] #A-CA-G
        epsilon_i_sd[37] = parameters[3]  #G-CA-A
        epsilon_i_sd[38] = parameters[42] #P-CA-G
        epsilon_i_sd[39] = parameters[40] #A-CT-G
        epsilon_i_sd[40] = parameters[4]  #G-CT-A
        epsilon_i_sd[41] = parameters[43] #P-CT-G
        epsilon_i_sd[42] = parameters[41] #A-CC-G
        epsilon_i_sd[43] = parameters[5]  #G-CC-A
        epsilon_i_sd[44] = parameters[44] #P-CC-G
        
        
        if model_id == 'seq_dep_finetuning_nearest_neighbour_both':
            return epsilon_si,rates_si,epsilon_i_sd
        elif model_id == 'seq_dep_finetuning_nearest_neighbour_both_constant_ei':
            epsilon_si[-20:] = 0.
            return epsilon_si,rates_si,epsilon_i_sd
        elif model_id == 'seq_dep_finetuning_nearest_neighbour_both_constant_ei_aba':
            epsilon_si_aba[-20:] = 0.
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
        else:
            return epsilon_si_aba,rates_si_aba,epsilon_i_sd
    
    else:
        print('Watch out! Non-existing model-ID..')
        return
    
    return epsilon,forward_rates

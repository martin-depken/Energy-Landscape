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
    
    if model_id == 'four_state_model_clv_engineered_cas2':
        if len(parameters)!=2:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = 2.242905038 #I
        epsilon[2] = parameters[0] #R
        epsilon[3:5] = 0. #barriers, set to zero, fit rates directly
        epsilon[5:25] = 6.
        
        epsilon[5]  -= 0.5 #1
        epsilon[6]  -= 3. #2
        epsilon[7]  -= 1. #3
        epsilon[8]  -= 0.1 #4
        epsilon[9]  -= 0.8 #5
        epsilon[10] -= 0.6 #6
        epsilon[11] -= 0.8 #7
        epsilon[12] -= 0.8 #8
        epsilon[13] += 0.2 #9
        epsilon[14] -= 2.5 #10
        epsilon[15] -= 1.5 #11
        epsilon[16] -= 2.5 #12
        epsilon[17] -= 1. #13
        epsilon[18] -= 1.4 #14
        epsilon[19] -= 1.3 #15
        epsilon[20] -= 0.8 #16
        epsilon[21] -= 1.7 #17
        epsilon[22] -= 3.0 #18
        epsilon[23] -= 2.5 #19
        epsilon[24] -= 4.3 #20
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1] = 10**-0.729266456785 #first rate
        forward_rates[2] = 10**parameters[1] #second rate
        forward_rates[3] = 10**1.9925159125999998 #cleavage rate
    
    elif model_id == 'four_state_model_clv_engineered_cas1':
        if len(parameters)!=2:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = 2.242905038 #I
        epsilon[2] = parameters[0] #R
        epsilon[3:5] = 0. #barriers, set to zero, fit rates directly
        epsilon[5:25] = 6.
        
        epsilon[5] -=0.5 #1
        epsilon[6] -=3. #2
        epsilon[7] -=1. #3
        epsilon[8] -= 0.1 #4
        epsilon[9] -=.5 #5
        epsilon[10] -=1.5 #6
        epsilon[11] -=2. #7
        epsilon[12] -=2. #8
        epsilon[13] += .5 #9
        epsilon[14] -= 1. #10
        epsilon[15] -= 1.5 #11
        epsilon[16] -= 1.3 #12
        epsilon[17] += 1. #13
        epsilon[18] += 2.5 #14
        epsilon[19] +=2.6 #15
        epsilon[20] += 2.5 #16
        epsilon[21] += 1.7 #17
        epsilon[22] += 1. #18
        epsilon[23] +=1.5 #19
        epsilon[24] -= 1.5 #20
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1] = 10**-0.729266456785 #first rate
        forward_rates[2] = 10**parameters[1] #second rate
        forward_rates[3] = 10**1.9925159125999998 #cleavage rate
    
    elif model_id == 'four_state_model_clv_engineered_cas':
        if len(parameters)!=2:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = 2.242905038 #I
        epsilon[2] = parameters[0] #R
        epsilon[3:5] = 0. #barriers, set to zero, fit rates directly
        epsilon[5:25] = [5.988477709800001, 4.48300119978, 5.63868077989, 6.19683729119, 5.1055995567900005, 6.146514016059999, 5.8447344994, 6.729104143, 7.49842519473, 5.69978855765, 5.5200466305499996, 4.5060459989699995, 5.2522936167, 5.10086442885, 4.53537972618, 3.73514366794, 2.50084758613, 2.5175495267400003, 3.06547635566, 2.50297376258]

        forward_rates[0] = 1000 #saturated
        forward_rates[1] = 10**-0.729266456785 #first rate
        forward_rates[2] = 10**parameters[1] #second rate
        forward_rates[3] = 10**1.9925159125999998 #cleavage rate
    
    
    
    elif model_id == 'four_state_model_clv_constant_ei':
        if len(parameters)!=8:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = parameters[0]
        epsilon[2] = parameters[1]
        epsilon[3:5] = parameters[2:4] #Barrier heights, relative to PAM, to I
        epsilon[5:21] = parameters[4] #Epsilon I
        epsilon[21:25] = parameters[5] #Epsilon I
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1:3] = 10**parameters[-2]
        forward_rates[3] = 10**parameters[-1]
        
    elif model_id == 'four_state_model_clv_rates_constant_ei':
        if len(parameters)!=6:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = parameters[0]
        epsilon[2] = parameters[1]
        epsilon[3:5] = 0. #Barrier heights, relative to PAM, to I
        epsilon[5:25] = parameters[2] #Epsilon I
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1] = 10**parameters[-3]
        forward_rates[2] = 10**parameters[-2]
        forward_rates[3] = 10**parameters[-1]
        
    elif model_id == 'four_state_model_on_rates_constant_ei':
        if len(parameters)!=7:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = parameters[0]
        epsilon[1] = parameters[1]
        epsilon[2] = parameters[2]
        epsilon[3:5] = 0. #Barrier heights, relative to PAM, to I
        epsilon[5:25] = parameters[3] #Epsilon I
        
        forward_rates[0] = 10**parameters[-3]
        forward_rates[1] = 10**parameters[-2]
        forward_rates[2] = 10**parameters[-1]
        forward_rates[3] = 0.
        
        
    elif model_id == 'engineered_clv_rates_constant_ei': #7_6_2019_sim_8_5, edited
        if len(parameters)!=3:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = 2.242905038
        epsilon[2] = parameters[0]
        epsilon[3:5] = 0. #Barrier heights, relative to PAM, to I
        epsilon[5:25] = 5.37 #Epsilon I
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1] = 10**-0.729266456785
        forward_rates[2] = 10**parameters[1]
        forward_rates[3] = 10**parameters[2]
        
    elif model_id == 'engineered_on_rates_constant_ei':#7_6_2019_sim_8_5, edited
        if len(parameters)!=2:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = 1.5641814477899998
        epsilon[1] = 2.242905038
        epsilon[2] = parameters[0]
        epsilon[3:5] = 0. #Barrier heights, relative to PAM, to I
        epsilon[5:25] = 5.37 #Epsilon I
        
        forward_rates[0] = 10**-2.39851897054
        forward_rates[1] = 10**-0.729266456785
        forward_rates[2] = 10**parameters[1]
        forward_rates[3] = 0.
        
    elif model_id == 'engineered_clv_rates_constant_ei_v2': #7_6_2019_sim_8_5, edited
        if len(parameters)!=4:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = parameters[0]
        epsilon[2] = parameters[1]
        epsilon[3:5] = 0. #Barrier heights, relative to PAM, to I
        epsilon[5:25] = 5.37 #Epsilon I
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1] = 10**-0.729266456785
        forward_rates[2] = 10**parameters[2]
        forward_rates[3] = 10**parameters[3]
        
    elif model_id == 'engineered_on_rates_constant_ei_v2':#7_6_2019_sim_8_5, edited
        if len(parameters)!=3:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = 1.5641814477899998
        epsilon[1] = parameters[0]
        epsilon[2] = parameters[1]
        epsilon[3:5] = 0. #Barrier heights, relative to PAM, to I
        epsilon[5:25] = 5.37 #Epsilon I
        
        forward_rates[0] = 10**-2.39851897054
        forward_rates[1] = 10**-0.729266456785
        forward_rates[2] = 10**parameters[2]
        forward_rates[3] = 0.
        
    elif model_id == 'four_state_model_on_rates_constant_ei_fixPAM':
        if len(parameters)!=6:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = 1.4
        epsilon[1] = parameters[0]
        epsilon[2] = parameters[1]
        epsilon[3:5] = 0. #Barrier heights, relative to PAM, to I
        epsilon[5:25] = parameters[2] #Epsilon I
        
        forward_rates[0] = 10**parameters[3]
        forward_rates[1] = 10**parameters[4]
        forward_rates[2] = 10**parameters[5]
        forward_rates[3] = 0.
    
    elif model_id == 'four_state_model_on_constant_ei':
        if len(parameters)!=9:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = parameters[0]
        epsilon[1] = parameters[1]
        epsilon[2] = parameters[2]
        epsilon[3:5] = parameters[3:5] #Barrier heights, relative to PAM, to I
        epsilon[5:21] = parameters[5] #Epsilon I
        epsilon[21:25] = parameters[6] #Epsilon I
        
        forward_rates[0] = 10**parameters[-2]
        forward_rates[1:3] = 10**parameters[-1]
        forward_rates[3] = 0.
        
    elif model_id == 'four_state_model_clv':
        if len(parameters)!=26:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = parameters[0]
        epsilon[2] = parameters[1]
        epsilon[3:5] = parameters[2:4] #Barrier heights, relative to PAM, to I
        epsilon[5:25] = parameters[4:24] #Epsilon I
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1:3] = 10**parameters[-2]
        forward_rates[3] = 10**parameters[-1]
    
    elif model_id == 'four_state_model_on':
        if len(parameters)!=27:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = parameters[0]
        epsilon[1] = parameters[1]
        epsilon[2] = parameters[2]
        epsilon[3:5] = parameters[3:5] #Barrier heights, relative to PAM, to I
        epsilon[5:25] = parameters[5:25] #Epsilon I
        
        forward_rates[0] = 10**parameters[-2]
        forward_rates[1:3] = 10**parameters[-1]
        forward_rates[3] = 0.
    
    elif model_id == 'four_state_model_clv_rates':
        if len(parameters)!=25:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = parameters[0]
        epsilon[2] = parameters[1]
        epsilon[3:5] = 0.
        epsilon[5:25] = parameters[2:22] #Epsilon I
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1] = 10**parameters[-3]
        forward_rates[2] = 10**parameters[-2]
        forward_rates[3] = 10**parameters[-1]
    
    elif model_id == 'four_state_model_on_rates':
        if len(parameters)!=26:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(5+20)
        forward_rates = np.ones(4)
        
        epsilon[0] = parameters[0]
        epsilon[1] = parameters[1]
        epsilon[2] = parameters[2]
        epsilon[3:5] = 0.
        epsilon[5:25] = parameters[3:23] #Epsilon I
        
        forward_rates[0] = 10**parameters[-3]
        forward_rates[1] = 10**parameters[-2]
        forward_rates[2] = 10**parameters[-1]
        forward_rates[3] = 0.
        
    elif model_id == 'six_state_model_clv':
        if len(parameters)!=26:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(7+20)
        forward_rates = np.ones(6)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = parameters[0]
        epsilon[2] = 0.
        epsilon[3] = parameters[1]
        epsilon[4] = 0.
        epsilon[5:7] = parameters[2:4] #Barrier heights, relative to PAM, to I
        epsilon[7:27] = parameters[4:24] #Epsilon I
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1:5] = 10**parameters[-2]
        forward_rates[5] = 10**parameters[-1]
        
        
    elif model_id == 'six_state_model_clv_constant_ei':
        if len(parameters)!=8:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(7+20)
        forward_rates = np.ones(6)
        
        epsilon[0] = -100 #saturated
        epsilon[1] = parameters[0]
        epsilon[2] = 0.
        epsilon[3] = parameters[1]
        epsilon[4] = 0.
        epsilon[5:7] = parameters[2:4] #Barrier heights, relative to PAM, to I
        epsilon[7:23] = parameters[4] #Epsilon I
        epsilon[23:27] = parameters[5] #Epsilon I
        
        forward_rates[0] = 1000 #saturated
        forward_rates[1:5] = 10**parameters[-2]
        forward_rates[5] = 10**parameters[-1]
    
    elif model_id == 'six_state_model_on':
        if len(parameters)!=27:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(7+20)
        forward_rates = np.ones(6)
        
        epsilon[0] = parameters[0]
        epsilon[1] = parameters[1]
        epsilon[2] = 0.
        epsilon[3] = parameters[2]
        epsilon[4] = 0.
        epsilon[5:7] = parameters[3:5] #Barrier heights, relative to PAM, to I
        epsilon[7:27] = parameters[5:25] #Epsilon I
        
        forward_rates[0] = 10**parameters[-2]
        forward_rates[1:5] = 10**parameters[-1]
        forward_rates[5] = 0.
        
    elif model_id == 'six_state_model_on_constant_ei':
        if len(parameters)!=9:
            print('Wrong number of parameters')
            return
        
        epsilon = np.zeros(7+20)
        forward_rates = np.ones(6)
        
        epsilon[0] = parameters[0]
        epsilon[1] = parameters[1]
        epsilon[2] = 0.
        epsilon[3] = parameters[2]
        epsilon[4] = 0.
        epsilon[5:7] = parameters[3:5] #Barrier heights, relative to PAM, to I
        epsilon[7:23] = parameters[5] #Epsilon I
        epsilon[23:27] = parameters[6] #Epsilon I
        
        forward_rates[0] = 10**parameters[-2]
        forward_rates[1:5] = 10**parameters[-1]
        forward_rates[5] = 0.
    
    else:
        print('Watch out! Non-existing model-ID..')
        return
    
    return epsilon,forward_rates

def combined_model(parameters,model_ID):
    
    model_ID_clv, model_ID_on = model_ID.split('+')
    
    if model_ID == 'six_state_model_clv+six_state_model_on':
        if len(parameters)!=28:
            print('Wrong number of parameters')
            return
        
        parameters_clv = np.append(parameters[1:25],parameters[26:28])
        parameters_on = np.array(parameters[0:27])
        
    elif model_ID == 'four_state_model_clv+four_state_model_on':
        if len(parameters)!=28:
            print('Wrong number of parameters')
            return
        
        parameters_clv = np.append(parameters[1:25],parameters[26:28])
        parameters_on = np.array(parameters[0:27])
        
    elif model_ID == 'four_state_model_clv_rates+four_state_model_on_rates':
        if len(parameters)!=27:
            print('Wrong number of parameters')
            return
        
        parameters_clv = np.append(parameters[1:23],parameters[24:27])
        parameters_on = np.array(parameters[0:26])
        
    elif model_ID == 'six_state_model_clv_constant_ei+six_state_model_on_constant_ei':
        if len(parameters)!=10:
            print('Wrong number of parameters')
            return
        
        parameters_clv = np.append(parameters[1:7],parameters[8:10])
        parameters_on = np.array(parameters[0:9])
        
    elif model_ID == 'four_state_model_clv_constant_ei+four_state_model_on_constant_ei':
        if len(parameters)!=10:
            print('Wrong number of parameters')
            return
        
        parameters_clv = np.append(parameters[1:7],parameters[8:10])
        parameters_on = np.array(parameters[0:9])
        
    elif model_ID == 'four_state_model_clv_rates_constant_ei+four_state_model_on_rates_constant_ei':
        if len(parameters)!=8:
            print('Wrong number of parameters')
            return
        
        parameters_clv = np.append(parameters[1:4],parameters[5:8])
        parameters_on = np.array(parameters[0:7])
        
    elif model_ID == 'four_state_model_clv_rates_constant_ei+four_state_model_on_rates_constant_ei_fixPAM':
        if len(parameters)!=7:
            print('Wrong number of parameters')
            return
        
        parameters_clv = np.append(parameters[0:3],parameters[4:7])
        parameters_on = np.array(parameters[0:6])
        
    elif model_ID == 'engineered_clv_rates_constant_ei+engineered_on_rates_constant_ei':
        if len(parameters)!=3:
            print('Wrong number of parameters')
            return
        
        parameters_clv = parameters
        parameters_on = parameters[0:2]
        
    elif model_ID == 'engineered_clv_rates_constant_ei_v2+engineered_on_rates_constant_ei_v2':
        if len(parameters)!=4:
            print('Wrong number of parameters')
            return
        
        parameters_clv = parameters
        parameters_on = parameters[0:3]
        
    else:
        print 'non-existing model-id'
        return
        
    return model_ID_clv, model_ID_on, parameters_clv, parameters_on
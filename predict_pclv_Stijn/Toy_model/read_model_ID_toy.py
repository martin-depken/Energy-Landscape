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
    if model_id == 'four_state_model_clv_constant_ei':
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
    else:
        print 'non-existing model-id'
        return
        
    return model_ID_clv, model_ID_on, parameters_clv, parameters_on
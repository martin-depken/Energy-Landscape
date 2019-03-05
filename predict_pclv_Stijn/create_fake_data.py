import numpy as np
import calculate_cleavage_rate as chi

def create_fake_data():
    #parameters = np.loadtxt('/home/svandersmagt/example-params-2.txt')
    parameters = np.loadtxt('example-params-2.txt')
    guide = 20
    model_ID = 'Clv_init_limit_Saturated_general_energies_v2'
    
    k = list()
    x = list()
    
    k.append(chi.calc_clv_rate_fast(parameters,model_ID,[],guide))
    x.append([])
    
    for i in range(1,guide+1):
        mismatch = [i]
        k.append(chi.calc_clv_rate_fast(parameters,model_ID,mismatch,guide))
        x.append(mismatch)
        
    for i in range(1,guide):
        for j in range(i+1,guide+1):
            mismatch = [i,j]
            k.append(chi.calc_clv_rate_fast(parameters,model_ID,mismatch,guide))
            x.append(mismatch)
            
    k = np.array(k)    
    error = k*0.05

    
    return x, k, error
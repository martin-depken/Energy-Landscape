from functions import *

#============================================
parameters = []
with open("example-params-2.txt",'r') as fl:
	for line in fl.readlines():
		parameters.append(float(line))
parameters = (parameters)
model_id = 'Clv_init_limit_general_energies_v2'


lookuptable_target = make_target_dict(parameters,model_id)
#print(lookuptable_target['A']['T'][19][0])
#print(lookuptable_target['A']['A'][1][0])
save_dict(lookuptable_target,"lookuptable_target")

lookuptable_PAM = make_PAM_dict(parameters,model_id)
save_dict(lookuptable_PAM,"lookuptable_PAM")
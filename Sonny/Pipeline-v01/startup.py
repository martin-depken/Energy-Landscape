""" This file should be executed on a usual computer."""

#============================================
__author__ = "Sonny Floyd de Jong"
__copyright__ = "Copyright 2019"
__license__ = "CC BY-NC-SA"
__version__ = "1.0.0"
__maintainer__ = "Depken lab"
__email__ = "s.f.dejong@student.tudelft.nl"
__status__ = "Production"
#============================================
from functions import *

#============================================
parameters = []

mainpath = open("mainpath.txt","r").readline()
if mainpath[-1] != "\\" :
	mainpath += '\\'
	
with open(mainpath+"example-params-2.txt",'r') as fl:
	for line in fl.readlines():
		parameters.append(float(line))
parameters = (parameters)
model_id = 'Clv_init_limit_general_energies_v2'


lookuptable_target = make_target_dict(parameters,model_id)
#print(lookuptable_target['A']['T'][1][0])
#print(lookuptable_target['A']['A'][1][0])
save_dict(lookuptable_target,mainpath+"lookuptable_target")

lookuptable_PAM = make_PAM_dict(parameters,model_id)
save_dict(lookuptable_PAM,mainpath+"lookuptable_PAM")
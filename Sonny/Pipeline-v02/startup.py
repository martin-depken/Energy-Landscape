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
import sys
from functions import *

#============================================
class CRISPR:
	def __init__(self, name):
		self.name = name
		if name in {'Cas9'}:
			self.guidelength = 20
			self._5primeseed_wrt_target = False
			self.PAM = True
			self.pamlength = 3
			self._5primePAM_wrt_target = False
			self.reversetarget = True
			self.complementtarget = True
		elif name in {'Cas12','Cas12a','Cpf1'}:
			self.guidelength = 20
			self._5primeseed_wrt_target = True
			self.PAM = True
			self.pamlength = 4
			self._5primePAM_wrt_target = True
			self.reversetarget = False
			self.complementtarget = False
			
def main(argv):
	Cas = CRISPR(argv[1])
	parameters = []

	try:
		mainpath = open("mainpath.txt","r").readline()
		if mainpath[-1] != "\\" :
			mainpath += '\\'
	except:
		mainpath = "/home/sfdejong/11_03_2019/"
	
	with open(mainpath+"parameters.txt",'r') as fl:
		for line in fl.readlines():
			parameters.append(float(line))
	parameters = (parameters)
	model_id = 'Clv_init_limit_general_energies_v2'
	
	with open(mainpath+"fit_22_3_2019_sim_13.txt",'r') as fl:
		lines = fl.read().splitlines()
		lastline = lines[-1].split("\t")
		print(type(lastline))
		parameters = [float(i) for i in lastline[:-3] ]
		
	parameters = tuple(parameters)
	print(parameters)
	model_id = 'general_energies_no_kPR'


	lookuptable_target = make_target_dict(parameters,model_id)
	print(lookuptable_target['A']['T'][1][0])
	print(lookuptable_target['A']['A'][1][0])
	save_dict(lookuptable_target,mainpath+"lookuptable_target")

	lookuptable_PAM = make_PAM_dict(parameters,model_id,Cas)
	save_dict(lookuptable_PAM,mainpath+"lookuptable_PAM")
	
if __name__ == "__main__":
	main(sys.argv)

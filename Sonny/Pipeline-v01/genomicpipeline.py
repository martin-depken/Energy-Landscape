""" This file contains a function that performs the backbone structure of the pipeline."""

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
#sys.path.append("/home/sfdejong/")
import hpc05notification
from functions import *
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
from time import time

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
	file = argv[1]
	Cas = CRISPR(argv[2])
	guide = "GGGUGGGGGGAGUUUGCUCC" #pVC297 VEGF Site#1 from 5' to 3'
	try:
		mainpath = open("mainpath.txt","r").readline()
		if mainpath[-1] != "\\" :
			mainpath += '\\'
	except:
		mainpath = "/home/sfdejong/11_03_2019/"

	t1 = time()
	lut_tar = read_dict(mainpath+"lookuptable_target")
	lut_pam = read_dict(mainpath+"lookuptable_PAM")
	startpos = 0 # PAM position at ChrY where it starts
	endpos = startpos+50*1187473
	
	if Cas.complementtarget == True: guide = str( Seq(guide).transcribe() )
	if Cas.reversetarget == True:    guide = guide[::-1]
	
	partition(file,mainpath,startpos,endpos,guide,lut_pam,lut_tar,Cas,includeNs=False)
	t2 = time()
	print("Elapsed time equals",t2-t1)
	try:
		hpc05notification.hpc05notification(t2-t1)
	except:
		print("Notification failed.")
	return
 
if __name__ == "__main__":
	main(sys.argv)

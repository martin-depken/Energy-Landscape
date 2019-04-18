""" This file contains a function that performs the backbone structure of the pipeline."""

#============================================
__author__ = "Sonny Floyd de Jong"
__copyright__ = "Copyright 2019"
__license__ = "CC BY-NC-SA"
__version__ = "2.0.0"
__maintainer__ = "Depken lab"
__email__ = "s.f.dejong@student.tudelft.nl"
__status__ = "Production"
#============================================

import sys
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
	Chr = argv[1]
	file = Chr+".fasta" # in the names of our filing system, this works
	Cas = CRISPR(argv[2])
	includePAMs = (argv[3]=="True")
	reverse = (argv[4]=="True") # false for positive strand
	guide = "GGGUGGGGGGAGUUUGCUCC" #pVC297 VEGF Site#1 from 5' to 3'
	
	try:
		mainpath = open("mainpath.txt","r").readline()
		if mainpath[-1] != "\\" :
			mainpath += '\\'
	except:
		mainpath = "/home/sfdejong/04_15_2019/"

	t1 = time()
	lut_tar = read_dict(mainpath+"lookuptable_target")
	lut_pam = read_dict(mainpath+"lookuptable_PAM")
	startpos = 0 # PAM position at ChrY where it starts
	
	if Cas.complementtarget == True: guide = str( Seq(guide).transcribe() )
	if Cas.reversetarget == True:    guide = guide[::-1]
	
	print(Chr)
	print("includePAMs = ",includePAMs)
	print("reverse = ",reverse)
	partition(file,mainpath,startpos,guide,lut_pam,lut_tar,Cas,reverse,Chr,False,includePAMs)
	t2 = time()
	print("Elapsed time equals",t2-t1)
	try:
		if (t2-t1) > 100:
			hpc05notification.hpc05notification(t2-t1,"final","comp")
	except:
		print("Notification failed.")
	return
 
if __name__ == "__main__":
	main(sys.argv)

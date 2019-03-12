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
from functions import *
from partition import *
import numpy as np
from Bio import Seq
from time import time


def main(argv):
	file = argv[1]
	mainpath = open("mainpath.txt","r").readline()
	if mainpath[-1] != "\\" :
		mainpath += '\\'
	
	t1 = time()
	lut_tar = read_dict(mainpath+"lookuptable_target")
	lut_pam = read_dict(mainpath+"lookuptable_PAM")
	startpos = 0 # PAM position at ChrY where it starts
	endpos = startpos+240
	guide = "GGGUGGGGGGAGUUUGCUCC" #pVC297 VEGF Site#1 from 5' to 3'
	partition(file,mainpath,startpos,endpos,guide,lut_pam,lut_tar,includeNs=False)
	t2 = time()
	print("Elapsed time equals",t2-t1)
	return
 
if __name__ == "__main__":
	main(sys.argv)

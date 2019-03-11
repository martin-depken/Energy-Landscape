import sys 
from functions import *
from partition import *
import numpy as np
from Bio.Seq import Seq


def main(argv):
	#parameter1 = argv[1]
	#print(parameter1, "algorithm started.")
	
	lut_tar = read_dict("lookuptable_target")
	lut_pam = read_dict("lookuptable_PAM")

	startpos = 10100 # PAM position at ChrY where it starts
	endpos = startpos+23
	if endpos<startpos+23:
		print("The end position is ill-defined.")
		
	guide = "GGGUGGGGGGAGUUUGCUCC" #pVC297 VEGF Site#1 from 5' to 3'
	guidelength = len(guide)
	file = "ChrY.fa"
	partition(file,startpos,endpos,guide,lut_pam,lut_tar,includeNs=False)
			
	return
 
if __name__ == "__main__":
	main(sys.argv)

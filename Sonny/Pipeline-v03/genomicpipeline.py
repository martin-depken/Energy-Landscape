""" This file contains a function that performs the backbone structure of the pipeline."""

#============================================
__author__ = "Sonny Floyd de Jong"
__copyright__ = "Copyright 2019"
__license__ = "CC BY-NC-SA"
__version__ = "3.0.0"
__maintainer__ = "Depken lab"
__email__ = "s.f.dejong@student.tudelft.nl"
__status__ = "Production"
#============================================

import sys
from functions import *
from Bio.Seq import Seq
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
	"""

	:param Chr: name of chromosome FASTA file without extension
	:param Cas: CRISPR/Cas type (see CRISPR class)
	:param includePAMs: indicates whether to include all PAMs or only NGG
	:param reverse: false for positive strand
	:param guide: RNA guide from 5' to 3' (irrespective of Cas type)
	:type Chr: str
	:type Cas: str
	:type includePAMs: bool
	:type reverse: bool

	:Example:
	>>>python genomicpipeline.py Chr01 Cas9 False False
	:return:
	"""

	# READING THE GIVEN ARGUMENTS
	Chr = argv[1]
	file = Chr+".fasta" # in the names of our filing system, this works
	Cas = CRISPR(argv[2]) # refers to the CRISPR class
	includePAMs = (argv[3]=="True") # indicates whether to include all PAMs or only NGG
	reverse = (argv[4]=="True") # false for positive strand
	if len(argv) == 6:
		guide = argv[5]
	else:
		guide = "GGGUGGGGGGAGUUUGCUCC" #pVC297 VEGF Site#1 from 5' to 3'

	# DEFINING THE MAIN PATH
	try:
		# on the cluster, this will not work
		mainpath = open("mainpath.txt","r").readline()
		if mainpath[-1] != "\\" :
			mainpath += '\\'
	except:
		mainpath = "/home/sfdejong/05_05_2019/"

	t1 = time()

	# READING LOOKUPTABLES
	lut_tar = read_dict(mainpath+"lookuptable_target")
	lut_pam = read_dict(mainpath+"lookuptable_PAM")

	startpos = 0 # at this position, the genome will be readPAM position at ChrY where it starts
	
	if Cas.complementtarget == True: guide = str( Seq(guide).transcribe() )
	if Cas.reversetarget == True:    guide = guide[::-1]
	
	print(Chr)
	print("includePAMs = ",includePAMs)
	print("reverse = ",reverse)
	print(guide)

	scanthegenome(file,mainpath,startpos,guide,lut_pam,lut_tar,Cas,reverse,Chr,False,includePAMs)

	t2 = time()
	print("Elapsed time equals",t2-t1)
	try:
		if (t2-t1) > 100 or True:
			hpc05notification.hpc05notification((t2-t1,str(Chr)+str(reverse)),"final","comp")
	except:
		print("Notification failed.")
	return
 
if __name__ == "__main__":
	main(sys.argv)

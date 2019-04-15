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
#sys.path.append("/home/sfdejong/")

from functions_anl import *
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
from time import time
from prediction_tools import *

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
	t1 = time()
	filename = argv[1]
	Cas = CRISPR(argv[2])
	#guide = "GGGUGGGGGGAGUUUGCUCC" #pVC297 VEGF Site#1 from 5' to 3'
	guide = "GGGTGGGGGGAGTTTGCTCC"
	PAM = "NGG"

	try:
		mainpath = open("mainpath.txt","r").readline()
		if mainpath[-1] != "\\" :
			mainpath += '\\'
	except:
		mainpath = "/home/sfdejong/03_21_2019/"
		filename = "/home/sfdejong/11_03_2019/chrY-B-data.hdf5"

	lut_tar = read_dict(mainpath+"lookuptable_target")
	lut_pam = read_dict(mainpath+"lookuptable_PAM")

	if filename[-5:] != ".hdf5":
		filename = filename+".hdf5"
	
	dsname = "ChrY"
	with h5py.File(filename, "r") as hdf:
		if not dsname in hdf.keys():
			print("Dataset (ChrY) did not exist.")
		#for longindex in range(hdf[dsname].shape[0]):
		score = np.array(hdf[dsname])
	
	steps = 100
	threshholds = np.linspace(np.min(score),np.max(score),steps)
	P = []
	N = []
	for sigma in threshholds:
		P.append(  np.sum(score >= sigma) )
		N.append(  np.sum(score < sigma) )

	genome_wide = pd.read_csv(mainpath+'genome_wide_precalc_models.csv')
	subset = genome_wide[(genome_wide['insertion_deletion'] == 'no') &
                         (genome_wide['truncated_guide'] == 'no') &
                         (genome_wide['canonical_PAM']=='yes') ]
						 
	targets  = subset[subset['guide']==guide+PAM]['target']
	guides  = subset[subset['guide']==guide+PAM]['guide']
	
	if Cas.complementtarget == True: guide = str( Seq(guide).transcribe() )
	if Cas.reversetarget == True:    guide = guide[::-1]	
	
	guide = str( Seq(guide).complement())
	truescores = []
	for seq in subset['target']:
		kclv = calculate_kclv(seq,guide,lut_pam,lut_tar,Cas,reverse=False)
		truescores.append( kclv[1])
	
	TP = []
	TN = []
	for sigma in threshholds:
		TP.append(  np.sum(truescores >= sigma) )
		TN.append(  np.sum(truescores < sigma) )
	print(truescores)
	print(TP,TN)
	
	O = pd.DataFrame()
	O['threshholds'] = threshholds
	O['P'] = P
	O['N'] = N
	O['TP'] = TP
	O['TN'] = TN
	O['FP'] = np.subtract(P , TP)
	O['FN'] = np.subtract(N , TN)
	O.to_csv("chrF-data.csv",index=True)
	
	return

if __name__ == "__main__":
	main(sys.argv)

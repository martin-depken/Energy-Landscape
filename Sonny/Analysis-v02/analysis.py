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
	calculategenomewide = True # false -> import the genome-wide data from a file

	try:
		mainpath = open("mainpath.txt","r").readline()
		if mainpath[-1] != "\\" :
			mainpath += '\\'
	except:
		mainpath = "/home/sfdejong/03_21_2019/"
		filename = mainpath + filename
	
	#lut_tar = read_dict("F:\\output\\lookuptable_target")
	#lut_pam = read_dict("F:\\output\\lookuptable_PAM")

	lut_tar = read_dict("lookuptable_target")
	lut_pam = read_dict("lookuptable_PAM")


	if filename[-5:] != ".hdf5":
		filename = filename+".hdf5"
	
	
	# ==============================
	# CALCULATE WHOLE-GENOME RESULTS
	# ==============================

	if calculategenomewide:
	
		steps = 1000
		prdP = np.zeros(steps) # predicted positives
		prdN = np.zeros(steps) # predicted negatives
		
		dsname = argv[1]
		#thresholds = np.linspace(0,2727688,steps)
		#thresholds = np.linspace(0,204741850000000,steps)
		#thresholds = np.linspace(5.2,1551786.5,steps)
		thresholds = np.linspace(24.387,7048387,steps) # minimal to maximal experimental pclv
		thresholds = np.linspace(24.387,12404500,steps)
		thresholds = np.linspace(2.5423722e-05,12404500,steps)
		
		maxfound = 0
		minfound = 1000
	
		with h5py.File(filename, "r") as hdf: # open HDF5 file
		
			for dsname in [key for key in hdf.keys()]: # loop over all datasets within this file
			
				print(dsname,[key for key in hdf.keys()])
				
				if int(np.ceil(hdf[dsname].shape[0]/1000000)) == 0:
					# if dataset is less than a million long
					
						score = np.array(hdf[dsname][0:None,1])
						maxfound = max(maxfound,np.max(score))
						minfound = min(minfound,np.min(score))

				else:
				
					for k in range(int(np.ceil(hdf[dsname].shape[0]/1000000))): # divide dataset into pieces of a million
					
						print("Now at: ",k)
					
						try:
							score = np.array(hdf[dsname][1000000*k:1000000*(k+1),1])
							maxfound = max(maxfound,np.max(score))
							minfound = min(minfound,np.min(score))
						except:
							score = np.array(hdf[dsname][1000000*k:None,1])
							maxfound = max(maxfound,np.max(score))
							minfound = min(minfound,np.min(score))
						
				for index in range(len(thresholds)):
					prdP[index] += np.sum(score >= thresholds[index])
					prdN[index] += np.sum(score <  thresholds[index])
		
		# print max and min from genome-wide file
		print("max:",maxfound)
		print("min:",minfound)
		
		print(prdP)
		print(prdN)
	
	else:
		
		O = pd.read_csv("analysis-results.csv")
		prdP = O['prdP']
		prdN = O['prdN']
		thresholds = O['thresholds']
	
	
	
	# ==============================
	# CALCULATE EXPERIMENTAL RESULTS
	# ==============================
	
	genome_wide = pd.read_csv(mainpath+'genome_wide_precalc_models.csv')
	subset = genome_wide[(genome_wide['Technique']== 'HTGTS') & (genome_wide['InDel'] == False) & ( (genome_wide['PAM'] =='AGG') | (genome_wide['PAM'] =='CGG') | (genome_wide['PAM'] =='GGG') | (genome_wide['PAM'] =='TGG')) &  (genome_wide['name']== 'VEGFA')]
	
	# select targets and guides
	targets = subset[subset['guide']==guide+PAM]['target']
	guides  = subset[subset['guide']==guide+PAM]['guide']
	
	if Cas.complementtarget == True: guide = str( Seq(guide).transcribe() )
	if Cas.reversetarget == True:    guide = guide[::-1]	
	
	guide = str( Seq(guide).complement())
	
	truescores = []
	print("Line 95")
	for seq in subset['target']:
		if 'N' in str(seq):
			print("Found an N!")
			continue
		kclv = calculate_kclv(seq,guide,lut_pam,lut_tar,Cas,reverse=False)
		truescores.append( kclv[0])
	print("Line 99")
	expP = np.zeros(len(thresholds))
	expN = np.zeros(len(thresholds))
	for index in range(len(thresholds)):
		expP[index] += np.sum(truescores >= thresholds[index])
		expN[index] += np.sum(truescores <  thresholds[index])
	
	O2 = pd.DataFrame()
	O2['thresholds'] = thresholds
	O2['prdP'] = prdP
	O2['prdN'] = prdN
	O2['expP'] = expP
	O2['expN'] = expN
	O2.to_csv("analysis-results.csv",index=True)
	
	return

if __name__ == "__main__":
	main(sys.argv)

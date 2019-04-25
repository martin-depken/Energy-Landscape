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
	t1 = time()
	filename = argv[1]
	Cas = CRISPR(argv[2])
	guide = "GGGUGGGGGGAGUUUGCUCC" #pVC297 VEGF Site#1 from 5' to 3'
	try:
		mainpath = open("mainpath.txt","r").readline()
		if mainpath[-1] != "\\" :
			mainpath += '\\'
	except:
		mainpath = "/home/sfdejong/03_21_2019/"
		filename = "/home/sfdejong/11_03_2019/chrY-B-data.hdf5"

	if filename[-5:] != ".hdf5":
		filename = filename+".hdf5"

	dsname = "ChrY"
	thresholds = np.array([200,4000,40000])
	thresholds = np.logspace(2,25,num=40,base=10,dtype='float')
	print(thresholds)
	positiveprediction = np.zeros([len(thresholds)])
	negativeprediction = np.zeros([len(thresholds)])
	with h5py.File(filename, "r") as hdf:
		if not dsname in hdf.keys():
			print("Dataset (ChrY) did not exist.")
		#for longindex in range(hdf[dsname].shape[0]):
		for longindex in range(20):
			if longindex%1000 == 0:
				stop = open(mainpath+"stop.txt","r").readline()
				if stop == "1":
					print("Stop!")
					break
			if longindex%round((hdf[dsname].shape[0]/5))==0:
				print("OK.")
				#hpc05notification.hpc05notification(longindex,"milestone","comp")
		
			tclv = hdf[dsname][longindex][1]
			
			if tclv < 0:
				negativeprediction += np.ones([len(thresholds)])
			else:
				test =  ( tclv < thresholds )*1
				#print(tclv,thresholds,test,1)
				positiveprediction += test
				negativeprediction += 1-test
			
	print(positiveprediction,negativeprediction)
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

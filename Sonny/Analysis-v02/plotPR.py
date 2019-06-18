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
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Raleway']})
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.style
import matplotlib
import matplotlib.cm


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
	#filename = argv[1]
	#Cas = CRISPR(argv[2])
	#guide = "GGGUGGGGGGAGUUUGCUCC" #pVC297 VEGF Site#1 from 5' to 3'
	guide = "GGGTGGGGGGAGTTTGCTCC"
	PAM = "NGG"
	
	# =============================

	O = pd.read_csv("analysis-results.csv")
	TP = O['expP']
	FN = O['expN'] #before drawing on paper
	TN = O['prdN'] - O['expN']
	FP = O['prdP'] - O['expP']	
		
	TPR  = TP.astype(float) / (TP + FN) # = recall = sensitivity
	FPR =  FP.astype(float) / (FP + TN)
	AUC_ROC = np.abs(np.trapz(y=TPR, x=FPR))

	#5) PR - curve:
	PRC =  TP.astype(float) / (TP + FP)
	# recall = TPR
	AUC_PR = np.abs(np.trapz(y=PRC[:-1], x=TPR[:-1]))	
	#print(P)
	
	
	O2 = pd.DataFrame()
	O2['thresholds'] = O['thresholds']
	O2['prdP'] = O['prdP']
	O2['prdN'] = O['prdN']
	O2['expP'] = O['expP']
	O2['expN'] = O['expN']
	O2['TP'] = TP
	O2['FN'] = FN
	O2['FP'] = FP
	O2['TN'] = TN
	O2['TPR']= TPR
	O2['FPR']= FPR
	O2['PRC']= PRC
	O2.to_csv("data-for-PR.csv",index=True)
	
	
	print("F1-score:", 2*PRC[50]*TPR[50]/(PRC[50]+TPR[50]))
	
	if True:
		plt.figure(figsize=(8, 6))
		ax = plt.subplot(111)    
		ax.spines["top"].set_visible(False)    
		ax.spines["bottom"].set_visible(False)    
		ax.spines["right"].set_visible(False)    
		ax.spines["left"].set_visible(False)
		ax.get_xaxis().tick_bottom()    
		ax.get_yaxis().tick_left() 

		hfont = {'fontname':'Raleway'}
		plt.plot(TPR,PRC,'x')
		plt.grid(False)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.title("PR curve",**hfont)
		plt.xlabel("TPR = RCL")
		plt.ylabel("PRC")
		plt.legend("x")
		ax.set_facecolor((0.9, 0.9, 0.9))
		ax.grid(color='w')
		#plt.savefig('example.pdf')
		plt.show()
		
		
		plt.figure(figsize=(8, 6))
		ax = plt.subplot(111)    
		ax.spines["top"].set_visible(False)    
		ax.spines["bottom"].set_visible(False)    
		ax.spines["right"].set_visible(False)    
		ax.spines["left"].set_visible(False)
		ax.get_xaxis().tick_bottom()    
		ax.get_yaxis().tick_left() 

		hfont = {'fontname':'Raleway'}
		plt.plot(FPR,TPR,'x')
		#plt.plot(FPR,TPR)
		plt.grid(False)
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		plt.title("ROC curve",**hfont)
		plt.xlabel("FPR")
		plt.ylabel("TPR")
		plt.legend("x")
		ax.set_facecolor((0.9, 0.9, 0.9))
		ax.grid(color='w')
		#plt.savefig('example.pdf')
		plt.show()	
		
	return

if __name__ == "__main__":
	main(sys.argv)

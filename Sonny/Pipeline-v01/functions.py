""" This file contains all the functions that are used elsewhere in the pipeline."""

#============================================
__author__ = "Sonny Floyd de Jong"
__copyright__ = "Copyright 2019"
__license__ = "CC BY-NC-SA"
__version__ = "1.0.0"
__maintainer__ = "Depken lab"
__email__ = "s.f.dejong@student.tudelft.nl"
__status__ = "Production"
#============================================

import pickle # for reading and saving dict
import numpy as np
from read_model_ID import *


#============================================

def save_dict(dict,name):
	"""Save a dictionary in the pickle format. The extension (.p) may be omitted.
	Keyword arguments:
	dict -- the (nested) dictionary)
	name -- the filename or file location (with or without .p extension)."""
	if name[-2:] == ".p":
		pickle.dump( dict, open(name, "wb"))
	else:
		pickle.dump( dict, open(name+".p", "wb"))
	return

def read_dict(name):
	"""Read a dictionary from a pickle format. The extension (.p) may be omitted.
	Keyword arguments:
	name -- the filename or file location (with or without .p extension)."""
	if name[-2:] == ".p":
		return pickle.load( open(name,"rb"))
	else:
		return pickle.load( open(name+".p","rb"))
		
def convert_param(parameters,model_id,guide,target):
	"""
	Uses the function unpack_parameters() to read a long list epsilons and forward rates. Next, it calculates the right epsilons for matching / mismatching nucleotides at each position and puts them in the following nested tuple:
	
	( (e_1,kf_1) , (e_2,kf_2) , ... , (e_20,kf_20) )
	
	Keyword arguments:
	parameters -- a list of parameters that unpack_parameters() can read
	model_id -- see unpack_parameters()
	guide -- one nucleotide from the guide RNA
	target -- one nucleotide from the target DNA
	
	Output:
	tuple in the form ( (e_1,kf_1) , (e_2,kf_2) , ... , (e_20,kf_20) )
	"""
	epsilon,forwardrates = unpack_parameters(parameters,model_id)
	matches = {'AA','CC','GG','UT'}
	if guide+target in matches:
		epsilons = -np.array(epsilon[1:21]) # epsilon_c (negative by contention)
	else:
		epsilons = -np.array(epsilon[1:21]) + np.array(epsilon[21:]) # epsilon_c + epsilon_i
	return tuple(zip(epsilons,forwardrates[2:]))
	
	
def generalise_dict(
	lut,
	equivalent=	[	('A','a','R','M','W','D','H','V','N'), 
					('C','c','S'), 
					('G','g'),
					['U','u','T','t','K','Y','B']
					]
					):
	"""
	Creates extra entries in the target lookup dict according to the lists stored in equivalent. This makes it possible to fastly extend the target lookup dict in case the FASTA format uses alternative nucleotides like R, M, W, D, H, V, S, K, Y, B etc.
	
	Keyword arguments:
	lut -- lookup dictionary for target nucleotides.
	equivalent -- list of tuples of keys. Each tuple starts with an existing key and continues with equivalent keys. Those will be appended to the dictionary.
	"""
	
	for equiv in equivalent:
		for eq in equivalent:
			letter = eq[0]
			if eq[0] == 'U': letter = 'T'
			for index in range(1,len(eq)):
				lut[ equiv[0] ][ eq[index] ] = lut[ equiv[0] ][ letter ]
		for index in range(1,len(equiv)):
			lut[ equiv[index] ] = lut[ equiv[0] ]
	return lut
	
def make_PAM_dict(
	parameters,
	model_id,
	equivalent=  [('A','a','N'),('C','c'),('G','g'),('U','u','T','t')]):
	"""
	Creates the PAM lookup dictionary using unpack_parameters().
	
	Keyword arguments:
	parameters -- list of parameters that unpack_parameters() can read
	model_id -- see unpack_parameters()
	equivalent -- list of tuples of nucleotides that are equivalent.
	
	Output:
	lut -- the PAM lookup dict with keys like 'NNG' that refer to tuples in the following form: (epsilon_PAM, forwardrate_PAM, forwardrate_solution).
	"""

	lut = {}
	epsilon,forwardrates = unpack_parameters(parameters,model_id)
	for letter1 in ('A','C','G','T','a','c','g','t'):
		for letter2 in ('A','C','G','T','a','c','g','t'):
			for letter3 in ('A','C','G','T','a','c','g','t'):
				lut[letter1+letter2+letter3] = (epsilon[0],forwardrates[1],forwardrates[0])
	return lut
	
	
def make_target_dict(parameters,model_id):
	"""
	Creates the target lookup dictionary using convert_param(). The output will be a nested dictionary with the following keys:
	
	epsilon_i,kf_i  = dict[guide][target][position=i]
	
	for instance,
	(3.93797145509, 235.7232903082969) =  dict['A']['T'][12]
	corresponds to the epsilon and forward rate for guide nucleotide 'A' and target nucleotide 'T' at position 12. Quite convenient, huh?
	
	Note that this function also generalises the dictionary before output.
	
	Keyword arguments:
	parameters --  list of parameters that unpack_parameters() can read
	model_id -- see unpack_parameters()
	"""
	lookuptable_target = {
	   #guide
		'A':	#target		#position in tuples of (epsilon,rate)
							# ( (e_1,kf_1) , (e_2,kf_2) , (e_3,kf_3) ... )
		
				{'A':		convert_param(parameters,model_id,'A','A'),
				 'C':		convert_param(parameters,model_id,'A','C'),
				 'G': 		convert_param(parameters,model_id,'A','G'),
				 'T': 		convert_param(parameters,model_id,'A','T')
				 },
		'C':
				{'A':		convert_param(parameters,model_id,'C','A'),
				 'C':		convert_param(parameters,model_id,'C','C'),
				 'G': 		convert_param(parameters,model_id,'C','G'),
				 'T': 		convert_param(parameters,model_id,'C','T')
				 },
		'G':
				{'A':		convert_param(parameters,model_id,'G','A'),
				 'C':		convert_param(parameters,model_id,'G','C'),
				 'G': 		convert_param(parameters,model_id,'G','G'),
				 'T': 		convert_param(parameters,model_id,'G','T')
				 },
		'U':
				{'A':		convert_param(parameters,model_id,'U','A'),
				 'C':		convert_param(parameters,model_id,'U','C'),
				 'G': 		convert_param(parameters,model_id,'U','G'),
				 'T': 		convert_param(parameters,model_id,'U','T')
				 },
		}
	return generalise_dict(lookuptable_target)
	

#============================================
# things below this line are not required for startup.py
#============================================

from Bio import SeqIO # for reading FASTA files
from Bio.Seq import Seq # reading DNA sequences
import h5py # handling HDF5 formats
import os

def readout_genome(name,mainpath,begin,end):
	"""
	This function will readout the genome using the 
	"""
	try:
		if name[-3:] == ".fa":
			gen = SeqIO.read(mainpath+"chromosomes/"+name, "fasta").seq
		elif name[-4:] == ".fna":
			gen = SeqIO.read(mainpath+"chromosomes/"+name, "fasta").seq
		else:
			print("Filename not recognised as FASTA format. Please add appropriate extension.")
	except ValueError:
		print("More than one record found in handle.")
		records = list(SeqIO.parse(mainpath+name, filetype))
		gen = records[0].seq
		for j in range(len(records)-1):
			gen += records[j+1].seq
	return gen[begin:end]

def single_target(target,guide,lut):
	"""
	Finds the energies and forward rates associated to a target sequence, guide sequence and lookup table.
	"""
	energy = np.zeros([20])
	rates = np.zeros([20])
	for index in range(20):
		energy[index],rates[index] = lut[ guide[index] ][ target[index] ][ 19-index ]
	return energy,rates
	
def store_tclv(position,startpos,tclv,filename,dsname="ChrY"):
	"""
	This function will store tclv in an hdf5 file.
	"""
	if filename[-5:] != ".hdf5":
		filename = filename+".hdf5"
	if not os.path.isfile(filename):
		print("File did not exist, so I created it.")
		with h5py.File(filename, "w-") as hdf:
			hdf.create_dataset(dsname, (1,2), maxshape=(None,2))
			hdf[dsname][0] = position+startpos,tclv
	else:
		with h5py.File(filename, "r+") as hdf:
			if not dsname in hdf.keys():
				print("Dataset (ChrY) did not exist, so I created it.")
				hdf.create_dataset(dsname, (0,2), maxshape=(None,2))
			dataindex = hdf[dsname].shape[0]
			hdf[dsname].resize((dataindex+1),axis=0)
			hdf[dsname][dataindex] = position+startpos,tclv
""" This file contains all the functions that are used elsewhere in the analysis."""

#============================================
__author__ = "Sonny Floyd de Jong"
__copyright__ = "Copyright 2019"
__license__ = "CC BY-NC-SA"
__version__ = "1.0.0"
__maintainer__ = "Depken lab"
__email__ = "s.f.dejong@student.tudelft.nl"
__status__ = "Production"
#============================================


import numpy as np
import h5py # handling HDF5 formats
import os
import hpc05notification
import pickle
from kinetic_model import *
from Bio.Seq import Seq

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
	
def read_HDF5(longindex,shortindex,filename,dsname="ChrY"):
	"""
	This function will store tclv in an hdf5 file.
	"""
	if filename[-5:] != ".hdf5":
		filename = filename+".hdf5"
	if not os.path.isfile(filename):
		print("File did not exist.")
		return

	with h5py.File(filename, "r") as hdf:
		if not dsname in hdf.keys():
			print("Dataset (ChrY) did not exist.")
			return
		return hdf[dsname][longindex][shortindex]
	return

def write_HDF5(longindex,shortindex,filename,value,dsname="ChrY"):
	"""
	This function will write a value in an hdf5 file.
	"""
	if filename[-5:] != ".hdf5":
		filename = filename+".hdf5"
	if not os.path.isfile(filename):
		print("File did not exist.")
		return
	if True:
		hdf = h5py.File(filename, "r+")
		if not dsname in hdf.keys():
			print("Dataset (ChrY) did not exist, so I created it.")
			hdf.create_dataset(dsname, (0,2), maxshape=(None,None))
		dataindex = hdf[dsname].shape[0]
		if longindex>dataindex:
			hdf[dsname].resize((dataindex+1),axis=0)
		print(longindex,shortindex,dsname,hdf[dsname].shape)
		print(hdf[dsname][0])
		print(hdf[dsname][longindex][shortindex])
		hdf[dsname][longindex][shortindex] = 20
		hdf.close()
	return

def read_dict(name):
	"""Read a dictionary from a pickle format. The extension (.p) may be omitted.
	Keyword arguments:
	name -- the filename or file location (with or without .p extension)."""
	if name[-2:] == ".p":
		return pickle.load( open(name,"rb"))
	else:
		return pickle.load( open(name+".p","rb"))
		
def single_target(target,guide,lut,Cas):
	"""
	Finds the energies and forward rates associated to a target sequence, guide sequence and lookup table.
	"""
	
		
	print(target,"t")
	print(guide,"g")
	
	energy = np.zeros([Cas.guidelength])
	rates = np.zeros([Cas.guidelength])
	for index in range(Cas.guidelength):
		energy[index],rates[index] = lut[ guide[index] ][ target[index] ][ index ]
	return energy,rates	
	
def calculate_kclv(sequence,guide,lut_pam,lut_tar,Cas,reverse=False):

	guidelength = Cas.guidelength
	pamlength   = Cas.pamlength
	_5prsd = Cas._5primeseed_wrt_target
	
	PAM = sequence[(not _5prsd)*(guidelength+pamlength - 1) :
			_5prsd*pamlength + (not _5prsd)*(guidelength-1):
			2*_5prsd-1]
	target = sequence[(guidelength-1)*(not _5prsd) +_5prsd*(pamlength): _5prsd*(guidelength+pamlength) or None :2*_5prsd-1]
	
	#calculate energies and forward rates
	pamenergy,pamrate,solrate = lut_pam[PAM]
	energy,forwardrates = single_target(target,guide,lut_tar,Cas)
	
	#concatenate and calculate backward rates
	energy = np.insert(energy,0,pamenergy)
	forwardrates = np.insert(forwardrates,0,[solrate,pamrate])

	backwardrates = get_backward_rates(energy, forwardrates,Cas)
	backwardrates = np.insert(backwardrates,0,[0])
	try:
		M = build_rate_matrix(forwardrates,backwardrates)
		tclv = mean_first_passage_time(M)
	except:
		M = build_rate_matrix(forwardrates,backwardrates)
		print("At position",position,"something went wrong with the matrix.\n",M)
		return "ERROR!"
	return max(1/(tclv*10**-9),0),tclv
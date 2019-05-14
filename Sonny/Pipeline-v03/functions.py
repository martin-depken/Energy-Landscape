""" This file contains all the functions that are used elsewhere in the pipeline."""

#============================================
__author__ = "Sonny Floyd de Jong"
__copyright__ = "Copyright 2019"
__license__ = "CC BY-NC-SA"
__version__ = "3.0.0"
__maintainer__ = "Depken lab"
__email__ = "s.f.dejong@student.tudelft.nl"
__status__ = "Production"
#============================================

import pickle # for reading and saving dict
from read_model_ID import *

from Bio import SeqIO # for reading FASTA files
import h5py # handling HDF5 formats
import os
from kinetic_model import *
import hpc05notification
from time import time

#============================================

def save_dict(dict,name):
	"""
	Save a dictionary in the pickle format. The extension (.p) may be omitted.

	:param dict: the (nested) dictionary)
	:param name: the filename or file location (with or without .p extension)
	:return:
	"""
	if name[-2:] == ".p":
		pickle.dump( dict, open(name, "wb"))
	else:
		pickle.dump( dict, open(name+".p", "wb"))
	return

def read_dict(name):
	"""
	Read a dictionary from a pickle format. The extension (.p) may be omitted.

	:param name: the filename or file location (with or without .p extension)
	:return:
	"""
	if name[-2:] == ".p":
		return pickle.load( open(name,"rb"))
	else:
		return pickle.load( open(name+".p","rb"))
		
def convert_param(parameters,model_id,guide,target):
	"""
	Uses the function unpack_parameters() to read a long list epsilons and forward rates. Next, it calculates the right epsilons for matching / mismatching nucleotides at each position and puts them in the following nested tuple:
	( (e_1,kf_1) , (e_2,kf_2) , ... , (e_20,kf_20) )

	:param parameters: a list of parameters that unpack_parameters() can read
	:param model_id: see unpack_parameters()
	:param guide: one nucleotide from the guide RNA
	:param target: one nucleotide from the target DNA
	:return: tuple in the form ( (e_1,kf_1) , (e_2,kf_2) , ... , (e_20,kf_20) )
	"""
	epsilon,forwardrates = unpack_parameters(parameters,model_id)
	forwardrates[-1] = 1000 # adjust the cleavage rate
	matches = {'AT','CG','GC','UA'}
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

	:param lut:  lookup dictionary for target nucleotides.
	:param equivalent: list of tuples of keys. Each tuple starts with an existing key and continues with equivalent keys. Those will be appended to the dictionary.
	:return:
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
	Cas,
	equivalent=  [('A','a','N'),('C','c'),('G','g'),('U','u','T','t')]):
	"""
	Creates the PAM lookup dictionary using unpack_parameters().

	:param parameters: list of parameters that unpack_parameters() can read
	:param model_id: see unpack_parameters()
	:param Cas:
	:param equivalent: list of tuples of nucleotides that are equivalent.
	:return: the PAM lookup dict with keys like 'NNG' that refer to tuples in the following form: (epsilon_PAM, forwardrate_PAM, forwardrate_solution).
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

	:param parameters: list of parameters that unpack_parameters() can read
	:param model_id: see unpack_parameters()
	:return:
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


def readout_genome(name,mainpath,begin,end=-1):
	"""
	This function will read out the DNA sequence from a FASTA file.

	:param name: name of FASTA file
	:param mainpath: location of FASTA file (eg: /home/sfdejong/chromosome/)
	:param begin:  first position from which the file will be read
	:param end:  last position from which the file will be read
	:return:
	"""
	try:
		if name[-3:] == ".fa":
			gen = SeqIO.read(mainpath+"chromosomes/"+name, "fasta").seq
		elif name[-4:] == ".fna":
			gen = SeqIO.read(mainpath+"chromosomes/"+name, "fasta").seq
		elif name[-6:] == ".fasta":
			gen = SeqIO.read(mainpath+"chromosomes/"+name, "fasta").seq
		else:
			raise ValueError("Filename not recognised as FASTA format. Please add appropriate extension.")
	except:
		raise ValueError("Something went wrong with reading out the FASTA file. Check readout_genome() and SeqIO.read()")

	if end != -1:
		# endpoint is defined
		return gen[begin:end]
	else:
		# end is not defined, so we take everything until the last nucleotide of the file
		return gen[begin:]

def get_energies_and_rates(target,guide,lut,Cas):
	"""
	Finds the energies and forward rates associated to a target sequence, guide sequence and lookup table.

	:param target: DNA target sequence from seed to end
	:param guide:  RNA (or DNA) guide sequence from seed to end
	:param lut:  lookuptable for target rates
	:param Cas:  object of CRISPR class
	:return: energies and forward rates
	"""
	energy = np.zeros([Cas.guidelength])
	rates = np.zeros([Cas.guidelength])
	for index in range(Cas.guidelength):
		energy[index],rates[index] = lut[ guide[index] ][ target[index] ][ index ]
	return energy,rates
	
def store_tclv(position,startpos,tclv,filename,pclv_min=0,dsname="chr"):
	"""
	This function will store tclv in an hdf5 file.

	:param position: position starting from position startpos
	:param startpos: first position of chromosome that has been read out
	:param tclv: average cleavage time
	:param filename: name of the HDF5 file
	:param pclv_min: value of pclv = 1/(tclv*10**-9) that acts as higher bound for too high tclv, where numerical errors occur
	:param dsname: prefix of the dataset (see file structure of HDF5)
	:return:
	"""
	if filename[-5:] != ".hdf5":
		filename = filename+".hdf5"

	# if the HDF5 file does not yet exist:
	if not os.path.isfile(filename):
		print("File did not exist, so I created it.")
		with h5py.File(filename, "w-") as hdf:
			hdf.create_dataset(dsname, (1,3), maxshape=(None,3))
			hdf[dsname][0] = position+startpos,max(1/(tclv*10**-9),pclv_min),tclv

	# if the HDF5 file already exists:
	else:
		with h5py.File(filename, "r+") as hdf:
			if not dsname in hdf.keys():
				print("Dataset (chr) did not exist, so I created it.")
				hdf.create_dataset(dsname, (0,3), maxshape=(None,3))
			dataindex = hdf[dsname].shape[0]
			hdf[dsname].resize((dataindex+1),axis=0)
			hdf[dsname][dataindex] = position+startpos,max(1/(tclv*10**-9),pclv_min),tclv
	return

def store_other_measures(position, startpos, guide, target, filename, mainpath, dsname="chr"):
	"""
	This function will store tclv in an hdf5 file.

	:param position: position starting from position startpos
	:param startpos: first position of chromosome that has been read out
	:param guide: guide RNA from seed to end
	:param target: target DNA from seed to end
	:param filename: name of the HDF5 file
	:param dsname: prefix of the dataset (see file structure of HDF5)
	:return:
	"""
	CFD,MIT = calculate_other_measures(guide,target, mainpath)

	if filename[-5:] != ".hdf5":
		filename = filename+".hdf5"

	# if the HDF5 file does not yet exist:
	if not os.path.isfile(filename):
		print("File did not exist, so I created it.")
		with h5py.File(filename, "w-") as hdf:
			hdf.create_dataset(dsname, (1,3), maxshape=(None,3))
			hdf[dsname][0] = position+startpos,CFD,MIT

	# if the HDF5 file already exists:
	else:
		with h5py.File(filename, "r+") as hdf:
			if not dsname in hdf.keys():
				print("Dataset (chr) did not exist, so I created it.")
				hdf.create_dataset(dsname, (0,3), maxshape=(None,3))
			dataindex = hdf[dsname].shape[0]
			hdf[dsname].resize((dataindex+1),axis=0)
			hdf[dsname][dataindex] = position+startpos,CFD,MIT
	return

def scanthegenome(file,mainpath,startpos,guide,lut_pam,lut_tar,Cas,reverse,Chr,includeNs=False,includePAMs=False):
	"""

	:param file: name of FASTA file with extension (eg: Chr01.fasta)
	:param mainpath: location of FASTA file (eg: /home/sfdejong/chromosome/)
	:param startpos: first position from which the file will be read
	:param guide: RNA sequence from seed to end, complementary to DNA in FASTA file
	:param lut_pam: lookuptable of PAM rates
	:param lut_tar: lookuptable of target rates
	:param Cas: object of CRISPR class
	:param reverse: false for positive strand
	:param Chr: name of chromosome FASTA file without extension
	:param includeNs:  indicates whether to include all undefined Ns or to skip them
	:param includePAMs: indicates whether to include all PAMs or only NGG
	:return:
	"""

	# DECLARE SOME VARIABLES
	guidelength = Cas.guidelength
	if guidelength != len(guide):
		raise ValueError("The given guide RNA sequence does not have a length that corresponds to the CRISPR/Cas protein.")
	pamlength   = Cas.pamlength
		
	# READ OUT GENOME FROM startpos UNTIL END (OR endpos)
	sequence = readout_genome(file,mainpath,startpos)
	endpos = len(sequence)

	# TAKE THE POSITIVE OR NEGATIVE STRAND OF THE CHROMOSOME
	if reverse:
		sequence = sequence.reverse_complement()
		filenamepart = "-"
	else:
		filenamepart = "+"

	# KEEP THE TIME
	lasttime = time()

	# SCAN THE PART OF THE GENOME THAT WE HAVE READ OUT
	for position in range(1,endpos-startpos-guidelength-pamlength):

		# EVERY 1000 POSITIONS, WE READ THE STOP FILE SO THAT WE CAN BREAK THE CALCULATIONS
		if position%1000 == 0:
			stop = open(mainpath+"stop.txt","r").readline()
			if stop == "1":
				print("Stop! We stopped the calculations, because the stop file (stop.txt) contained a 1.")
				break

		# SEND AN EMAIL ON EVERY QUARTER OF THE CHROMOSOME
		if position%(int(endpos/4)) == 0:
			if (time()-lasttime)>1200:
				try:
					hpc05notification.hpc05notification((str(position),str(endpos),str(Chr)+str(filenamepart)),"milestone","comp")
				except:
					print("Notification failed.")

		# READOUT PAM AND target IN THE CORRECT ORIENTATION, DEPENDING ON THE Cas OBJECT
		_5prsd = Cas._5primeseed_wrt_target
		PAM = sequence[
			position + (not _5prsd)*(guidelength+pamlength - 1) :
			position + _5prsd*pamlength + (not _5prsd)*(guidelength-1):
			2*_5prsd-1]
		target = sequence[
			position + _5prsd*(pamlength) + (not _5prsd)*(guidelength - 1) :
			(position + _5prsd*(guidelength+pamlength) -1*(not _5prsd) ) :
			2*_5prsd-1]

		# these PAMs will be admitted
		admissiblePAMS = {'GGA','GGC','GGG','GGT'}

		# WE WILL ONLY PERFORM CALCULATIONS IF PAM+target DOES NOT CONTAIN OTHER NUCLEOTIDES THAN ACGT (and N)
		if str(PAM+target).translate(None,'ACGTN') != "":
			print("Found an unknown nucleotide in ",str(PAM+target)," at position ",position)
			continue

		# CALCULATE AND SAVE THE AVERAGE CLEAVAGE TIME
		if ((not 'N' in PAM+target) or includeNs) and ((PAM in admissiblePAMS) or includePAMs):

			# print(PAM+"|"+target+str(position))

			# calculate energies and forward rates
			pamenergy,pamrate,solrate = lut_pam[PAM]
			energy,forwardrates = get_energies_and_rates(target,guide,lut_tar,Cas)
			
			# concatenate and calculate backward rates
			energy = np.insert(energy,0,pamenergy)
			forwardrates = np.insert(forwardrates,0,[solrate,pamrate])
			backwardrates = get_backward_rates(energy, forwardrates,Cas)
			backwardrates = np.insert(backwardrates,0,[0])
			
			# build rate matrix and calculate tclv			
			try:
				M = build_rate_matrix(forwardrates,backwardrates)
				tclv = mean_first_passage_time(M)
			except:
				M = build_rate_matrix(forwardrates,backwardrates)
				print("At position",position,"something went wrong with the matrix.\n",M)
				continue

			# save tclv or pclv into an HDF5 file
			store_tclv(position,startpos,tclv,mainpath+"output/"+Chr+filenamepart+".hdf5",0,Chr+filenamepart)
			try:
				store_other_measures(position, startpos, guide, target, mainpath+"output/alt-"+Chr+filenamepart+".hdf5", mainpath, Chr+filenamepart)
			except:
				print("XAt position",position,"something went wrong with the other measures.")

	return
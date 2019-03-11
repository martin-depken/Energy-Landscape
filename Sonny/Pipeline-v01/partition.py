import numpy as np
from functions import *
from kinetic_model import *

def partition(file,startpos,endpos,guide,lut_pam,lut_tar,includeNs=False):
	guidelength = len(guide)

	if endpos<startpos+3+guidelength:
		print("The end position is ill-defined.")
		
	sequence = readout_genome("chromosomes/"+file,startpos,endpos)


	for position in range(endpos-startpos-guidelength-2):
		PAM = sequence[position:position+3]
		target = sequence[position+3: position+3+20]
		if (not 'N' in PAM+target) or includeNs:
			target = "CCCACCCCCCTCAAACGAGG"
			
			#calculate energies and forward rates
			pamenergy,pamrate,solrate = lut_pam[PAM]
			energy,forwardrates = single_target(target,guide,lut_tar)
			
			#concatenate and calculate backward rates
			energy = np.insert(energy,0,pamenergy)
			forwardrates = np.insert(forwardrates,0,[solrate,pamrate])
			backwardrates = get_backward_rates(energy, forwardrates)
			backwardrates = np.insert(backwardrates,0,[0])
			
			# build rate matrix and calculate tclv			
			try:
				M = build_rate_matrix(forwardrates,backwardrates)
				tclv = mean_first_passage_time(M)
				print(position+startpos,PAM+"|"+target,tclv,position)
			except:
				M = build_rate_matrix(forwardrates,backwardrates)
				print("At position",position,"something went wrong with the matrix.\n",M)
				continue
			
			store_tclv(position,startpos,tclv,"dataset001.hdf5")
			
	return
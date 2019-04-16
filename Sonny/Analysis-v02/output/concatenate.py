""" This file contains a function that puts all hdf5 chromosome files in one file."""

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
import os
import numpy as np
import h5py

def main(argv):
	try:
		mainpath = open("mainpath.txt","r").readline()
		if mainpath[-1] != "\\" :
			mainpath += '\\'
	except:
		mainpath = "/home/sfdejong/04_15_2019/output/"

	joinedfile = mainpath+"joined.hdf5"
	with h5py.File(joinedfile, "w") as mainhdf:
		for file in os.listdir(mainpath):
			if file.endswith(".hdf5") and not file.endswith("joined.hdf5"):
				print("->", mainpath+file,"now concatenated.")
				with h5py.File(file, "r") as hdf:
					dsname = [key for key in hdf.keys()][0]
					print(dsname)
					hdf.copy(dsname,mainhdf)
					print("Done?")
	return

if __name__ == "__main__":
	main(sys.argv)

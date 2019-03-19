# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:55:58 2019

@author: Diewertje
"""
# Function to load the final parameters that got out of the fit

def load_simm_anneal(filename, Nparams):
    f=open(filename)
    lines=f.read().splitlines()
    last_line=lines[-1]
    last_line = last_line.split()
    parameters = list(map(float,last_line[:-2]))
    return parameters
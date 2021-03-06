# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 12:28:32 2019

@author: Diewertje
To make lists with the MMposition, ABA and uncertainty with every measurement. 
The input is the full data file, output is the lists with the information we need for the simulations.
For delta ABA data
"""
import numpy as np 
import pandas as pd


# path = '../Data_ABA_Finkelsteinlab/' 
# filename= 'cas9-target-e-replicate-1-delta-abas_Canonical_OT-r_0-20.csv'
def Prepare_Cdata(path, filename):
    # function berouz to get this csv file!
    
    data=pd.read_csv(path+filename)
    
    Grouped = data.groupby('Mutation Positions').agg(lambda x: list(x))
    Grouped.reset_index(inplace=True)
    
    for i in range(len(Grouped)):
        s=Grouped['Mutation Positions'][i]
        Grouped['Mutation Positions'][i] = np.array(s.split('|')).astype(int)
    
    MMpos=Grouped['Mutation Positions'].tolist()
    ABA=Grouped['Delta ABA (kBT)'].tolist()
    Uncertainty=Grouped['Uncertainty'].tolist()
    
    return MMpos, ABA, Uncertainty


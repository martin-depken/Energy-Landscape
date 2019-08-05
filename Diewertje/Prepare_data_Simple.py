# -*- coding: utf-8 -*-
"""
Created on Fri Feb 29 15:43:32 2019

@author: Diewertje
"""
import numpy as np 
import pandas as pd


# path = '../Data_ABA_Finkelsteinlab/' 
# filename= 'cas9-target-e-replicate-1-delta-abas_Canonical_OT-r_0-20.csv'
def Prepare_Cdata(path, filename):
    # function berouz to get this csv file!
    
    data=pd.read_csv(path+filename)
    
    Grouped = data.groupby('Mutation Positions').agg(lambda x: list(x))
    # In doing this grouping, we exclude the ontarget values, this did not matter for the deltaABA, but now it actually does
    # The ontarget measurements should have been included in the fits. 
    Grouped.reset_index(inplace=True)
    
    for i in range(len(Grouped)):
        s=Grouped['Mutation Positions'][i]
        Grouped['Mutation Positions'][i] = np.array(s.split('|')).astype(int)
    
    MMpos=Grouped['Mutation Positions'].tolist()
    ABA=Grouped['ABA'].tolist()
    Uncertainty=Grouped['error'].tolist()
    
    return MMpos, ABA, Uncertainty


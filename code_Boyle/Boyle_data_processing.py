import numpy as np
import pandas as pd


''''
Preprosses the data files from Boyle et al. to be used for Simmulated Annealing

Behrouz Eslami & Misha Klein    Depken lab 
'''


def prepare_multiprocessing(replica='1', path='../Data_Boyle/',
                            use_on_rate=True, use_off_rate=True, use_occupancy=True,
                            use_blocks_only=False):
    '''
    Prepares the data in such a format that it is usable for multiprocessing and Simmulated Annealing.

    To (hopefully) optimally make use of the multiprocessing during the fitting we return:
    xdata = [ [1,1] , [1,2] , [],.....,[20,20]  ] with every element specifying the two (or one) mismatch position(s)

    ydata = [ (occ_data, on_data,off_data) , (occ_data,on_data,off_data), ....]
    with every element of ydata containing a set of three datasets
    corresponding to the mismatch configuration shown of the coreesponding element of xdata

    yerr = [ (occ_err, on_err,off_err) , (occ_err,on_err,off_err), ....]
    with every element of ydata containing a set of three datasets
    corresponding to the mismatch configuration shown of the coreesponding element of xdata

    :param replica: 1 or 2
    :param path: path to data files
    :return:
    '''

    # 1) Read the datafiles and preposses them into a Pandas dataframe:
    Combined_Data = read(replica, path, consecutive_mm=use_blocks_only)

    # 2) Get the three datasets (if we choose to not fit any of the datasets, supply empty lists in stead.
    # The Chi2 calculation then knows how to handle it) 
    if use_occupancy:
        occ = np.array(Combined_Data.occ)
        occ_error = np.array(Combined_Data.occ_error)
    else:
        occ = []
        occ_error = []
        for i in range(len(Combined_Data.occ)):
            occ.append([])
            occ_error.append([])
        occ = np.array(occ)
        occ_error = np.array(occ_error)

    if use_on_rate:
        on = np.array(Combined_Data.on_slope)
        on_error = np.array(Combined_Data.on_error)
    else:
        on = []
        on_error = []
        for i in range(len(Combined_Data.on_slope)):
            on.append([])
            on_error.append([])
        on = np.array(on)
        on_error = np.array(on_error)

    if use_off_rate:
        off = np.array(Combined_Data.off_slope)
        off_error = np.array(Combined_Data.off_error)
    else:
        off = []
        off_error = []
        for i in range(len(Combined_Data.off_slope)):
            off.append([])
            off_error.append([])
        off = np.array(off)
        off_error = np.array(off_error)

    # 3) Get ydata in the format wanted
    ydata = []
    for i in range(len(occ)):
        ydata.append((occ[i], on[i], off[i]))

    # 4) Get xdata (Mismatch locations)
    xdata = np.array(Combined_Data.MM_pos)

    # 5) Get yerr in the format wanted
    yerr = []
    for i in range(len(occ)):
        yerr.append((occ_error[i], on_error[i], off_error[i]))
    return xdata, ydata, yerr

def weights_averages(replica='1', path='../Data_Boyle/'):
    '''
    Weights used for different datasets (occupancy, on-rate and off-rate)
    within the Chi-Square calculation during fitting

    Here use the averages of each dataset as the weights to get all terms in Chi-Squared
    of same order of magnitude (use inverse square for ChiSquare).

    :param replica: ID of '1' or '2'
    :param path: data files Boyle (not collected)
    :return: vector with the weights to be fed into the Simulated Annealing code
    '''
    DataTable = read(replica,path)
    OccWeight = (DataTable.occ.apply(np.mean).mean() )**(-2)
    OnWeight = (DataTable.on_slope.apply(np.mean).dropna().mean() )**(-2)
    OffWeight= (DataTable.off_slope.apply(np.mean).dropna().mean())**(-2)
    return np.array([OccWeight,OnWeight,OffWeight])


def read(replica='1',path='../Data_Boyle/', consecutive_mm=False):
    '''
    Load the datafiles with the original data from Boyle et al. and aggregate them into a single dataframe
    :param replica: 1 or 2
    :param path: path to data files
    :return: pandas dataframe
    '''
    occ = pd.read_csv(path + 'occupancy_rep' + str(replica) + '_processed_with_errors.txt')
    kon = pd.read_csv(path + 'second_fit_data.summarized.on.1nM.rep' + str(replica) + '.txt', delimiter='\t')
    koff = pd.read_csv(path + 'second_fit_data.summarized.off.10nM.rep' + str(replica) + '.txt', delimiter='\t')
    kon = kon[kon['nmut'] < 3][['mutations', 'slope', 'se']]
    koff = koff[koff['nmut'] < 3][['mutations', 'slope', 'se']]
    occ = occ[['mutations', 'Ratio', 'rel_err_Ratio']]
    kon.rename(columns={'slope': 'on_slope'}, inplace=True)
    koff.rename(columns={'slope': 'off_slope'}, inplace=True)
    occ.rename(columns={'Ratio': 'occ'}, inplace=True)
    kon.rename(columns={'se': 'on_error'}, inplace=True)
    koff.rename(columns={'se': 'off_error'}, inplace=True)
    occ.rename(columns={'rel_err_Ratio': 'occ_error'}, inplace=True)
    koff['off_slope'] *= -1
    Full_data = kon.copy()
    Full_data = Full_data.merge(koff, how='outer', on='mutations')
    Full_data = Full_data.merge(occ, how='outer', on='mutations')
    Full_data['MM_pos'] = Full_data['mutations'].apply(get_pos)
    Full_data['PAM_mut'] = Full_data['MM_pos'].apply(find_PAM_mutations)
    Full_data['mm_ID'] = Full_data['MM_pos'].apply(string_get_pos)
    no_PAM = Full_data[Full_data.PAM_mut == False]
    Combined_Data = no_PAM[['mm_ID', 'occ', 'on_slope', 'off_slope', 'on_error', 'off_error', 'occ_error']]
    Combined_Data = Combined_Data.groupby('mm_ID').agg(lambda x: list(x))
    Combined_Data = Combined_Data.reset_index()
    Combined_Data['MM_pos'] = Combined_Data['mm_ID'].apply(get_pos_again)
    Combined_Data = Combined_Data[['MM_pos', 'occ', 'on_slope', 'off_slope', 'on_error', 'off_error', 'occ_error']]
    Combined_Data['occ'] = Combined_Data['occ'].apply(remove_NaN)
    Combined_Data['on_slope'] = Combined_Data['on_slope'].apply(remove_NaN)
    Combined_Data['off_slope'] = Combined_Data['off_slope'].apply(remove_NaN)
    Combined_Data['occ_error'] = Combined_Data['occ_error'].apply(remove_NaN)
    Combined_Data['on_error'] = Combined_Data['on_error'].apply(remove_NaN)
    Combined_Data['off_error'] = Combined_Data['off_error'].apply(remove_NaN)


    if consecutive_mm:
        Combined_Data = Combined_Data[Combined_Data['MM_pos'].apply(find_consecutive)]
    return(Combined_Data)


'''
Helper functions 
'''
def find_consecutive(mm_pos):
    '''
    Only keep those mismatch configurations with consecutive mismatches, single mismatches or the on-target
    :param mm_pos:
    :return:
    '''
    return(np.abs(np.diff(mm_pos)) == 1).all() or (len(mm_pos)<2)






def get_pos(s):
    if s == 'WT':
        return([])
    s2 = s.split(':')
    return map(lambda x: 20-int(x[:-1]), s2)

def find_PAM_mutations(x):
    return (np.array(x)<=0).any()

def string_get_pos(x):
    s = ''
    for pos in x:
        s += str(pos) + 'X:'
    if len(x)==0:
        s='WT:'
    return s[:-1]

def get_pos_again(s):
    if s == 'WT':
        return([])
    s2 = s.split(':')
    return np.array(map(lambda x: int(x[:-1]), s2)).astype('int')

def remove_NaN(x):
    return np.array(x)[np.where(np.isnan(x)==False)[0]]
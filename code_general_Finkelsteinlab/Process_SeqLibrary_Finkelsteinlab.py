import numpy as np
import pandas as pd
import Bio
from Bio.pairwise2 import format_alignment
import copy

'''
Main functions
'''

def process_Finkelstein_Library(data, on_target, seq_colname, output_colnames, Cas='Cas9', Canonical_PAM = True,
                               Mut_type = ['OT','r'], Mut_min = 0, Mut_max = 2, out_file_name = 'data',
                                out_path = './', nofilter = False, save_data = True):
    '''

    This function processes the sequence library used in Finkelstein lab and produces and saves a data frame with all
    relevant information.

    :param data: raw data (data frame)
    :param on_target: the sequence of the intended target (non-target strand, 27 bp long)
    :param seq_colname: the name of the column in the data frame which contains the off-targets
    :param output_colnames: The names of the columns in the data frame that we want to keep (list of strings)
    :param Cas: Type of Cas protein: ['Cas9', 'Cas12a']
    :param Canonical_PAM: only keep off-targets with canonical PAM
    :param Mut_type: The type of mutations to be kept in the data frame (['OT', 'r', 'i', 'd'])
    :param Mut_min: Minimum number of mutations to be kept in the data frame.
    :param Mut_max: Maximum number of mutations to be kept in the data frame.
    :param out_file_name: The name of the file in which the data frame is saved
    :param out_path: the Path where the data frame is saved
    :param nofilter: If True, ignores the filtering criteria and returns the full data set
    :param save_data:  if True, saves the data frame
    :return:
    The processed data frame. It contains the off-target sequence column as well as the columns given in
    output_colnames from the original data set. Moreover, it contains the following columns:
    'On Target': Whether the off-target is the same as on-target. NOTE: this only includes the hybrid part
                 (PAM excluded). The PAM could be canonical or not.
    'PAM': PAM sequence.
    'Canonical': True if PAM is canonical. False otherwise.
    'Length difference': Difference in length between off-target and on-target.
    'Alignment': A string which shows the alignment between the two sequence (use print to visualize nicely).
                 NOTE: PAM is excluded. In the case of Cas9, the sequence is inverted. Therefore for both Cas9
                 and Cas12a, the left-mpst base is the most closest to PAM.
                 On-target sequence shown on top, off-target sequence in the bottom.
                 Empty string if undetermined.
    'Alignment (All)': All possible allignments discovered by Biopython.
    'Mutation Type' 'r': replacement, 'i': insertion, 'd': deletiom. Empty string if undetermined.
    'Mutation Count': Number of mutations (PAM excluded). NA if undetermined.
    'Mutation Positions': A string identifying the position of mutations, seperated by '|':
                          'pos1|pos2|pos3...'
                          Empty string if undetermined.
                          NOTE: PAM is excluded, positions starts with 1, and the smaller p[ositions are
                          closer to PAM. This is directly usable in our codes.
    'Mutation ID': A string showing the type, the sequence and the position of each mutation, for example:
                   'r:10:C|d:13:T|i:17:T|....'.
                   Empty string if undetermined.

    '''

    new_data = data.copy()
    cols_tokeep = [seq_colname]+output_colnames
    new_data = new_data[cols_tokeep]
    new_data['On Target'] = new_data[seq_colname].apply(lambda S: find_ontarget(S, on_target, Cas))
    new_data['PAM'] = new_data[seq_colname].apply(lambda S: separate_PAM(S, Cas)[0])
    new_data['Canonical'] = new_data[seq_colname].apply(lambda S: separate_PAM(S, Cas)[2])
    new_data['Length difference'] = new_data[seq_colname].apply(lambda S: find_length_diff(S, on_target, Cas))
    new_data['Alignment_raw'] = new_data[seq_colname].apply(lambda S: Allign(S, on_target, Cas))
    new_data['Alignment_selected'] = new_data.apply(Clean_allignment, axis=1)
    new_data['Alignment'] = new_data['Alignment_selected'].apply(lambda A: format_alignment(*A[0]).split('S')[0] if len(A) > 0 else '')
    new_data['Alignment (All)'] = new_data['Alignment_raw'].apply(lambda A: ('\n').join(map(lambda a: format_alignment(*a).split('S')[0], A)) if len(A)>0 else '')
    new_data['Mutation ID'] = new_data.apply(lambda x: Make_Mutation_ID(x,seq_colname,on_target, Cas), axis=1)
    new_data['Mutation Type'] = new_data['Mutation ID'].apply(lambda x: '|'.join(list(set(map(lambda y: y[0], x.split('|'))))) if ((x!='')&(x!='OT')) else '')
    new_data['Mutation Count'] = new_data['Mutation ID'].apply(lambda x: len(x.split('|')) if ((x!='')&(x!='OT')) else np.NaN)
    new_data.loc[new_data['Mutation ID']=='OT','On Target'] = True
    new_data.loc[new_data['On Target'],'Mutation Count'] = 0
    new_data.loc[new_data['On Target'],'Mutation Type'] = 'OT'
    new_data['Mutation Positions'] = new_data['Mutation ID'].apply(lambda x: '|'.join(map(lambda y: y.split(':')[1],x.split('|'))) if ((x!='')&(x!='OT')) else '')
    new_data.drop(['Alignment_raw','Alignment_selected'], inplace=True, axis=1)
    select_Pam = new_data['Canonical']== Canonical_PAM
    select_muttype = new_data['Mutation Type'].apply(lambda x: x in Mut_type)
    select_mutcount = new_data['Mutation Count'].apply(lambda x: (x>=Mut_min)&(x<=Mut_max))
    select = select_Pam&select_muttype&select_mutcount
    if not nofilter:
        new_data = new_data[select]
    PAM_type = 'nonCanonical'
    if Canonical_PAM:
        PAM_type = 'Canonical'
    Mut_str = '-'.join(Mut_type)
    output_filename = out_path+out_file_name+'_'+PAM_type+'_'+Mut_str+'_'+str(Mut_min)+'-'+str(Mut_max)+'.csv'
    if nofilter:
        output_filename = out_path+out_file_name+'_'+'full'+'.csv'
    if save_data:
        new_data.to_csv(output_filename, index=False)
    return new_data


'''
Helper functions 
'''

def separate_PAM(S, Cas='Cas9'):
    if Cas == 'Cas9':
        PAM_len = 3
        PAM = S[-PAM_len:]
        s = S[0:-PAM_len]
        s = s[::-1]
        s = s[0:-4]
        canonical = False
        if PAM[1:] == 'GG':
            canonical = True
    if Cas == 'Cas12a':
        PAM_len = 4
        PAM = S[:PAM_len]
        s = S[PAM_len:]
        s = s[0:-3]
        canonical = False
        if (PAM[:3] == 'TTT') & (PAM[3] != 'T'):
            canonical = True
    return PAM, s, canonical


def find_length_diff(S, on_target, Cas='Cas9'):
    _, s, canonical = separate_PAM(S, Cas)
    _, t, _ = separate_PAM(on_target, Cas)
    return (len(s) - len(t))


def find_ontarget(S, on_target, Cas='Cas9'):
    _, s, _ = separate_PAM(S, Cas)
    _, t, _ = separate_PAM(on_target, Cas)
    return s == t


def Allign(S, on_target, Cas='Cas9'):
    _, s, _ = separate_PAM(S, Cas)
    _, t, _ = separate_PAM(on_target, Cas)
    A = Bio.pairwise2.align.globalxx(t, s)
    if len(s) == len(t):
        MM_num = np.sum(np.array(list(t)) == np.array(list(s)))
        a = (t, s, MM_num, 0, len(t))
        if a not in A:
            A.append(a)

    return A


def Clean_allignment(x):
    A = copy.deepcopy(x['Alignment_raw'])
    Length_diff = x['Length difference']

    if Length_diff == 0:
        function_to_filter = lambda a: not (('-' in a[0]) or ('-' in a[1]))
        A = filter(function_to_filter, A)
        return A

    if Length_diff > 0:
        function_to_filter = lambda a: not ((a[0].count('-') != np.abs(Length_diff)) or ('-' in a[1]))
        A = filter(function_to_filter, A)
        id_for_sort = lambda a: '|'.join(
            map(lambda x: str(x), list(np.arange(1, len(a[0]) + 1)[np.array(list(a[0])) != np.array(list(a[1]))])))
        if len(A) > 1:
            A.sort(key=id_for_sort, reverse=True)
            del A[1:]
        return (A)

    if Length_diff < 0:
        function_to_filter = lambda a: not (('-' in a[0]) or (a[1].count('-') != np.abs(Length_diff)))
        A = filter(function_to_filter, A)
        id_for_sort = lambda a: '|'.join(
            map(lambda x: str(x), list(np.arange(1, len(a[0]) + 1)[np.array(list(a[0])) != np.array(list(a[1]))])))
        if len(A) > 1:
            A.sort(key=id_for_sort, reverse=True)
            del A[1:]
        return (A)


def Make_Mutation_ID(x, seq_colname, on_target, Cas='Cas9'):
    _, t, _ = separate_PAM(on_target, Cas)
    _, s, _ = separate_PAM(x[seq_colname], Cas)

    if x['On Target']:
        return 'OT'

    if len(x['Alignment_selected']) == 0:
        return ''

    Length_diff = int(x['Length difference'])
    a = x['Alignment_selected'][0]

    if Length_diff > 0:
        ta = np.array(list(a[0]))
        sa = np.array(list(a[1]))
        All_positions = np.arange(1, len(sa) + 1)
        Mut_positions_on_s = All_positions[ta != sa]
        Mut_Seqs = sa[[ta != sa]]
        Mut_types = np.array(['i'] * len(Mut_Seqs))
        Mut_Seqs_on_t = ta[[ta != sa]]
        Mut_types[Mut_Seqs_on_t != '-'] = 'r'
        offset = np.zeros(len(Mut_types), dtype=int)
        for n in range(len(offset)):
            offset[n] = np.sum(Mut_types[0:n] == 'i')
        # For insertion it gives the position of a base on the target before which there is an insertion.
        Mut_positions = Mut_positions_on_s - offset
        ID_list = []
        for Mut_type, Mut_pos, Mut_Seq in zip(Mut_types, Mut_positions, Mut_Seqs):
            ID_list.append(':'.join([Mut_type, str(Mut_pos), Mut_Seq]))
        ID_list = filter(lambda x: 'i:21' not in x[:-2], ID_list)  # Insertion after position 20 does not count.
        ID = '|'.join(ID_list)
        if ((ID == '') & (s[:-Length_diff] == t)):
            ID = 'OT'

    if Length_diff < 0:
        ta = np.array(list(a[0]))
        sa = np.array(list(a[1]))
        All_positions = np.arange(1, len(ta) + 1)
        Mut_positions = All_positions[ta != sa]
        Mut_Seqs = ta[[ta != sa]]
        Mut_types = np.array(['d'] * len(Mut_Seqs))
        Mut_Seqs_on_s = sa[[ta != sa]]
        Mut_types[Mut_Seqs_on_s != '-'] = 'r'
        Mut_Seqs[Mut_Seqs_on_s != '-'] = Mut_Seqs_on_s[Mut_Seqs_on_s != '-']
        ID_list = []
        for Mut_type, Mut_pos, Mut_Seq in zip(Mut_types, Mut_positions, Mut_Seqs):
            ID_list.append(':'.join([Mut_type, str(Mut_pos), Mut_Seq]))
        ID = '|'.join(ID_list)

    if Length_diff == 0:
        ta = np.array(list(a[0]))
        sa = np.array(list(a[1]))
        All_positions = np.arange(1, len(ta) + 1)
        Mut_positions = All_positions[ta != sa]
        Mut_Seqs = sa[[ta != sa]]
        Mut_types = np.array(['r'] * len(Mut_Seqs))
        ID_list = []
        for Mut_type, Mut_pos, Mut_Seq in zip(Mut_types, Mut_positions, Mut_Seqs):
            ID_list.append(':'.join([Mut_type, str(Mut_pos), Mut_Seq]))
        ID = '|'.join(ID_list)

    return ID

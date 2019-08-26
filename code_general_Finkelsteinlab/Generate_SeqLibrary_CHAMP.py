import pandas as pd
import numpy as np
import pickle
import Process_SeqLibrary_Finkelsteinlab as Finkellib
from scipy.optimize import curve_fit
import matplotlib.pylab as plt

def Generate_SeqLibrary(filename,
                        concentrations,
                        nmbr_boot,
                        on_target_seq,
                        non_specific_seq='AAGGCCGAATTCTCACCGGCCCCAAGGTATTCAAG',
                        mode='ref_sequence',
                        Cas='Cas9'):

    raw_CHAMP = unpickle_datadict(filename)

    df_CHAMP = collect_all_clusters(raw_CHAMP, concentrations, on_target_seq,Cas)

    del raw_CHAMP
    SeqLibrary = build_dataframe_for_fitting(df_CHAMP,
                                             nmbr_boot,
                                             concentrations,
                                            non_specific_seq,
                                            mode)
    return SeqLibrary



def unpickle_datadict(filename):
    data = pickle.load(open(filename))
    datadict = dict(data)
    return datadict


def collect_all_clusters(CHAMP_data, concentrations, target_seq, Cas):
    '''
    returns a dataframe with every row being a single cluster.
    Hence, every sequence is repeated as there are multiple clusters with the same sequence


    Can be used later to get the median binding curve
    '''
    Seq = []
    nmbr_concentration_points = len(CHAMP_data[CHAMP_data.keys()[0]][0])
    if nmbr_concentration_points != len(concentrations):
        print "wrong number of concentration points. Length of concentration vector does not match number of measurements"
        return

    c = [[] for i in range(nmbr_concentration_points)]

    for sequence in CHAMP_data.keys():
        for cluster, curve in enumerate(CHAMP_data[sequence]):
            Seq.append(sequence)
            for i in range(len(concentrations)):
                c[i].append(curve[i])
    df = pd.DataFrame()
    df['sequence'] = Seq
    for j in range(len(concentrations)):
        df[str(concentrations[j]) + 'nM'] = c[j]

    # find on-target, canonical PAM , etc.
    df['On Target'] = df['sequence'].apply(lambda S: Finkellib.find_ontarget(S, target_seq, Cas))
    df['PAM'] = df['sequence'].apply(lambda S: Finkellib.separate_PAM(S, Cas)[0])
    df['Canonical'] = df['sequence'].apply(lambda S: Finkellib.separate_PAM(S, Cas)[2])
    return df



def build_dataframe_for_fitting(full_data, nmbr_boot, concentrations,
                                non_specific_seq='AAGGCCGAATTCTCACCGGCCCCAAGGTATTCAAG',
                                mode='ref_sequence'):
       # remove the clusters for which there is not a complete binding curve available:
    full_data.dropna(axis=0, inplace=True)

    # --- substract background ----
    if mode == 'ref_sequence':
        minus_Imin = subtract_background_ref_sequence(full_data, non_specific_seq, concentrations)
    elif mode == 'lowest_concentration':
        minus_Imin = subtract_background_low_concentration(full_data, concentrations)

    # --- saturation level (to convert into bound fraction) ----
    Imax = find_saturation(minus_Imin, concentrations)

    # --- for every sequence: use bootstrapping to find avg. ABA and errorbar ABA ----
    sequences = []
    ABA = []
    error = []

    all_seqs = minus_Imin.drop(['On Target', 'PAM', 'Canonical'], axis=1).copy()

    for name, group in all_seqs.groupby('sequence'):
        # --- perform bootstrapping on ABA values for this sequence  ----
        avg, std = bootstrap_ABA(group, nmbr_boot, Imax, concentrations)

        # --- store results into processed dataframe ----
        sequences.append(name)
        ABA.append(avg)
        error.append(std)

    # --- stroe results into processed dataframe ----
    processed_data = pd.DataFrame()
    processed_data['sequence'] = sequences
    processed_data['ABA'] = ABA
    processed_data['error'] = error
    return processed_data


def get_binding_curves(df, concentrations, non_specific_seq='AAGGCCGAATTCTCACCGGCCCCAAGGTATTCAAG',mode='ref_sequence'):
    df.dropna(axis=0, inplace=True)
    # --- substract background ----
    if mode == 'ref_sequence':
        minus_Imin = subtract_background_ref_sequence(df, non_specific_seq, concentrations)
    elif mode == 'lowest_concentration':
        minus_Imin = subtract_background_low_concentration(df, concentrations)
    Imax = find_saturation(minus_Imin, concentrations)

    binding_curves = minus_Imin.drop(['On Target', 'PAM', 'Canonical'], axis=1).copy()

    median_binding_curves = binding_curves.groupby('sequence').median()

    # median_binding_curves.set_index('sequence', inplace=True)
    median_binding_curves /= Imax
    return median_binding_curves


#####################################################################################
def bootstrap_ABA(grouped, nmbr_runs, Imax, concentrations):
    bstrp_ABA = []
    if nmbr_runs>0:
        for run in range(nmbr_runs):
            # --- bootstrap sample cluster (per given sequence) -----
            bstrp_sample = grouped.sample(n=len(grouped), replace=True)

            # --- get our new 'median' binding curve ----
            bstrp_med = bstrp_sample.median()
            bstrp_binding_curve = row_to_list(bstrp_med, concentrations) / Imax
            # --- calculate ABA ---
            try:
                bstrp_ABA.append(np.log(curve_fit(Hill_eq, concentrations, bstrp_binding_curve, maxfev=1000)[0][0]))
            except RuntimeError:
                bstrp_ABA.append(np.nan)
    else:
        bstrp_med = grouped.median()
        bstrp_binding_curve = row_to_list(bstrp_med, concentrations) / Imax
        try:
            bstrp_ABA.append(np.log(curve_fit(Hill_eq, concentrations[concentrations>0],
                                              bstrp_binding_curve[concentrations>0], maxfev=1000)[0][0]))
        except RuntimeError:
            bstrp_ABA.append(np.nan)
    bstrp_ABA = np.array(bstrp_ABA)
    return np.nanmean(bstrp_ABA), np.nanstd(bstrp_ABA)


def subtract_background_low_concentration(df, concentrations):
    new_df = df.copy()
    background = df[str(concentrations[0]) + 'nM'].copy()
    for c in concentrations:
        new_df[str(c) + 'nM'] = new_df[str(c) + 'nM'] - background
    return new_df


def subtract_background_ref_sequence(df, non_specific_seq, concentrations):
    new_df = df.copy()
    subset_control_seq = new_df[new_df.sequence == non_specific_seq]
    background = subset_control_seq[str(concentrations[0]) + 'nM'].median()
    for c in concentrations:
        new_df[str(c) + 'nM'] = new_df[str(c) + 'nM'] - background
    return new_df


def find_saturation(df, concentrations):
    OnTarget = df[df['On Target'] & df['Canonical']].copy()
    saturation_value = OnTarget[str(concentrations[-1]) + 'nM'].median()
    return saturation_value


def row_to_list(row, concentrations):
    binding_curve = []
    for c in concentrations:
        binding_curve.append(row[str(c) + 'nM'])
    return np.array(binding_curve)
    
def Hill_eq(C, Kd):
    return (1.0 + Kd / C) ** (-1)


def get_ABA(df, concentrations):
    df['ABA'] = df['binding_curves'].apply(lambda x: np.log(curve_fit(Hill_eq, concentrations, x)[0][0]))
    return df
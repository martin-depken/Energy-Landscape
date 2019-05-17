import numpy as np

def get_backward_rates(energies, forwardrates, Cas):
	# kb(n) = kf(n-1) * exp(+epsilon(n))
    backwardrates = np.zeros(Cas.guidelength)
    backwardrates = forwardrates[:-1] * np.exp(energies)
    return backwardrates

def build_rate_matrix(forward_rates, backward_rates):
    diagonal1 = - (forward_rates + backward_rates)
    diagonal2 = backward_rates[1:]
    diagonal3 = forward_rates[:-1]
    rate_matrix = np.diag(diagonal1, k=0) + np.diag(diagonal2, k=1) + np.diag(diagonal3, k=-1)
    return rate_matrix

def mean_first_passage_time(M):
    guidelength = len(M)
    everything_unbound = np.array([1.0] + [0.0] * (guidelength-1))
    return np.sum(np.linalg.solve(-M,everything_unbound))

###
# BELOW THIS, WE HAVE OTHER MEASURES
###

def seq2num(gen):
    # Set constants
    # Please notice we assign Adenine to any nucleotide for which its nature is undetermined
    dict = {'A' : 0, 'a' : 0,
            'C' : 1, 'c': 1,
            'G' : 2, 'g' : 2,
            'T' : 3, 't' : 3,
            'N' : 0, 'n' : 0,
            'U' : 3, 'u' : 3}
    Nt = len(gen)  # Length (bp) of target genome

    # Initialize vector
    num_gen = np.zeros(Nt)
    # Convert sequence of nucleotides into a sequence of numbers
    for j in range(Nt):
        num_gen[j] = dict[gen[j]]

    return num_gen.astype(int)

import pandas as pd
def calculate_other_measures(guide,target,mainpath):

    # CALCULATE CFD
    cfd_table = pd.read_csv(mainpath+
        'STable 19 FractionActive_dlfc_lookup.csv')

    # CONVERT SEQUENCES TO NUMBERS
    guide = seq2num(guide)
    target = seq2num(target)

    DNA = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    RNA = {0: 'A', 1: 'C', 2: 'G', 3: 'U'}

    mm_pos = np.where(guide.astype(int)[:20] != target.astype(int)[:20])[0]
    # table shows the target strand (complementary), we search for homology
    target_strand = 3 - target
    score = 1.

    for MMpos in mm_pos:
        MMtype = 'r' + RNA[int(guide[MMpos])] + ':d' + DNA[int(target_strand[MMpos])]
        CFD = cfd_table[(cfd_table['Mismatch Type'] == MMtype) & (cfd_table['Position'] == MMpos + 1)]['Percent-Active'].iloc[0]
        score *= CFD
    print(score)
    CFD = score



    M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804,
         0.685, 0.583]
    mm_pos = np.where(guide.astype(int)[:20] != target.astype(int)[:20])[0]
    # mm_pos = np.where(guide != target)[0]
    n_mm = len(mm_pos)
    if n_mm > 0:
        dtot = 0
        for pos1 in range(0, n_mm):
            for pos2 in range(pos1 + 1, n_mm):
                dtot += abs(mm_pos[pos1] - mm_pos[pos2])
        if n_mm > 1:
            davg = dtot / (n_mm * (n_mm - 1) / 2.)
        else:
            davg = 0.

        factor_distance = 1 / ((19 - davg) / 19 * 4 + 1) * 1. / n_mm ** 2

        factor_positions = 1
        for pos in mm_pos:
            factor_positions *= (1 - M[pos])

        score = factor_distance * factor_positions
    else:
        score = 1.0
    MIT=score
    return CFD,MIT
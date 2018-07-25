# Boyle, Andreasson, Chircus, et al.
# Written by Evan Boyle
# eaboyle@stanford.edu

import pandas as pd 
import numpy as np
from heapq import *
from random import random
import argparse
import sys
from multiprocessing import Pool
from collections import defaultdict
from scipy.optimize import minimize
from joblib import Parallel, delayed
from copy import deepcopy

#lambda1_target = "GACGCATAAAGATGAGACGC"
#egfp_site1_target = "GGGCACGGGCAGCTTGCCGG"
#egfp_site2_target = "GATGCCGTTCTTCTGCTTGT" 

temperature = 298.15
R = 8.3145 

from math import exp,log10
from string import maketrans
from scipy.spatial.distance import euclidean
comp = maketrans("ACGT", "TGCA")


parser = argparse.ArgumentParser(description='Simulate sgRNA invasion')

parser.add_argument('-i', '--input_sequence', type=str, required=True,
                    help='DNA target sequence (nontemplate strand)')
parser.add_argument('-t', '--nn_table', type=str, required=True, 
                   help='path to file of thermodynamic parameters')
parser.add_argument('-f', '--association_rate', type=float, default=0.000425,
                   help='rate from state 0 to state 1')
parser.add_argument('-g', '--halt_fraction', type=float, required=True,
                   help='halts once provided fraction bound is reached in wild type')
parser.add_argument('-v','--transversion_penalty', type=float, default=11000,
                   help='energy penalty for Tv mismatch')
parser.add_argument('-s','--transition_penalty', type=float, default=7500,
                   help='energy penalty for Ts mismatch')
parser.add_argument('-b','--pam_binding_energy', type=float, default=-5000,
                   help='delta G for PAM binding')
parser.add_argument('-a','--energy_adjustment', type=float, default=-26000,
                   help='energy adjustment for protein purposes')
parser.add_argument('-p','--adjustment_position', type=int, default=12,
                   help='energy adjustment position for protein purposes')
parser.add_argument('-o', '--out_prefix', type=str, required=True,
                   help='output prefix')
parser.add_argument('-c', '--cluster_size', type=int, default=1000,
                   help='Number of DNAs to simulate')
parser.add_argument('-k', '--length_of_fixed_time', type=float,
                    help='Run simulation for fixed pseudotime')
parser.add_argument('-l', '--length_of_max_time', type=float, default=40000.,
                   help='Max pseudotime to run')
parser.add_argument('-n', '--n_simulations', type=int, default=1,
                   help='Number of clusters')
parser.add_argument('-d', '--dissociation', action='store_true',
                   help='Simulate dissociation rather than association')
parser.add_argument('-r', '--reference_table', type=str, required=True,
                   help='Rates for calculating correlation')
parser.add_argument('-w', '--wt_slope', type=float, default=1.,
                   help='Value to normalize slope')
parser.add_argument('-m', '--time_multiplier', type=float, default=10,
                    help='Time relative to run (in association time units) before simulating dissociation')
parser.add_argument('-u', '--n_wt_simulations', type=int, default=5,
                   help='Number of wt clusters for time estimation')
parser.add_argument('-j', '--n_jobs', type=int, default=3,
                   help='Number of cores to run in parallel')
parser.add_argument('-x', '--exit', action="store_true",
                   help='Exit after loading args for interactive purposes')

args = parser.parse_args()

def rc(seq):
    return seq[::-1].translate(comp)

nn_data = pd.read_table(args.nn_table)
wt_target = rc(args.input_sequence)
slope_data = pd.read_table(args.reference_table)

class dna:
    def __init__(self, wt_target, ref_dGs, mut_states, mut_bases, args, state=0):
        self.mut_states = mut_states
        self.mut_bases = mut_bases
        self.dGs = get_mutant_dGs(\
            wt_target, \
            ref_dGs, \
            args.transition_penalty, \
            args.transversion_penalty, \
            mut_states, \
            dict(zip(mut_states, [rc(mut_base) for mut_base in mut_bases])))

        self.state = state
        self.v_fs = np.array([get_v_f_dG(i, self.dGs) for i in range(0,21)])
        self.v_fs[0] = args.association_rate
        self.v_rs = np.array([get_v_r_dG(i, self.dGs) for i in range(0,21)])
        self.next_state = -1
        self.time_step = -1
        self.delta = 0

    def move_and_prepare(self):
        self.delta = 1. if self.state == 0 and self.next_state == 1 else -1. if self.state == 1 and self.next_state == 0 else 0.
        self.state = self.next_state
        self.prepare_move()
        return(self.time_step)

    def prepare_move(self):
        self.time_step = -log10(random()) / (self.v_fs[self.state] + self.v_rs[self.state])
        self.next_state = (self.state + 1) if random() < (self.v_fs[self.state] / (self.v_fs[self.state] + self.v_rs[self.state])) else (self.state - 1)    

        return(self.time_step)


# state, mut positions will be 1-based
def get_delta_G(state, mut_positions, base_lookup, pam_binding_energy):
    RNA_init = 0
    DNA_init = 0
    RNA_base_stack = 0
    DNA_base_stack = 0
    PAM_dG = pam_binding_energy if state > 0 else 0

    for i in range(1,state):
        if i in mut_positions:
            RNA_data = nn_data.loc[(nn_data.MISMATCH_BASE == "None") & \
            (nn_data.NA_1 == "RNA") & \
            (nn_data.NA_2 == "DNA") & \
            (nn_data.BASE_A == rc(wt_target[i])) & \
            (nn_data.BASE_B == rc(base_lookup[i])),]
        elif ((i + 1) in mut_positions):
            RNA_data = nn_data.loc[(nn_data.MISMATCH_BASE == "None") & \
            (nn_data.NA_1 == "RNA") & \
            (nn_data.NA_2 == "DNA") & \
            (nn_data.BASE_A == rc(base_lookup[i + 1])) & \
            (nn_data.BASE_B == rc(wt_target[i - 1])),]
        else:
            RNA_data = nn_data.loc[(nn_data.MISMATCH_BASE == "None") & \
            (nn_data.NA_1 == "RNA") & \
            (nn_data.NA_2 == "DNA") & \
            (nn_data.BASE_A == rc(wt_target[i])) & \
            (nn_data.BASE_B == rc(wt_target[i - 1])),]
       
        RNA_base_stack = RNA_base_stack + float(RNA_data.dH - temperature * RNA_data.dS)
    
    if state > 1:
        RNA_start_data = nn_data.loc[(nn_data.NA_1 == "RNA") & \
            (nn_data.NA_2 == "DNA") & \
            (nn_data.BASE_A == "ATGC") & \
            (nn_data.BASE_B == "init"),]
        RNA_init = RNA_init + float(RNA_start_data.dH - temperature * RNA_start_data.dS)
        RNA_end_data = nn_data.loc[(nn_data.NA_1 == "RNA") & \
            (nn_data.NA_2 == "DNA") & \
            (nn_data.BASE_A == "ATGC") & \
            (nn_data.BASE_B == "init"),]
        RNA_init = RNA_init + float(RNA_end_data.dH - temperature * RNA_end_data.dS)

    
    for i in range(1 if state == 1 else state + 1, 20): # state 1 is PAM binding, no melting
        if i in mut_positions:
            DNA_data = nn_data.loc[(nn_data.MISMATCH_BASE == "None") & \
                (nn_data.NA_1 == "DNA") & \
                (nn_data.NA_2 == "DNA") & \
                (nn_data.BASE_A == base_lookup[i]) & \
                (nn_data.BASE_B == wt_target[i]),]
        elif ((i + 1) in mut_positions):
            DNA_data = nn_data.loc[(nn_data.MISMATCH_BASE == "None") & \
                (nn_data.NA_1 == "DNA") & \
                (nn_data.NA_2 == "DNA") & \
                (nn_data.BASE_A == wt_target[i - 1]) & \
                (nn_data.BASE_B == base_lookup[i + 1]),]
        else:
            DNA_data = nn_data.loc[(nn_data.MISMATCH_BASE == "None") & \
                (nn_data.NA_1 == "DNA") & \
                (nn_data.NA_2 == "DNA") & \
                (nn_data.BASE_A == wt_target[i - 1]) & \
                (nn_data.BASE_B == wt_target[i]),]

        DNA_base_stack = DNA_base_stack + float(DNA_data.dH - temperature * DNA_data.dS)

    if wt_target[19] == "A" or wt_target[19] == "T":
        DNA_end_data = nn_data.loc[(nn_data.NA_1 == "DNA") & \
            (nn_data.NA_2 == "DNA") & \
            (nn_data.BASE_A == "AT") & \
            (nn_data.BASE_B == "init"),]
    else:
        DNA_end_data = nn_data.loc[(nn_data.NA_1 == "DNA") & \
            (nn_data.NA_2 == "DNA") & \
            (nn_data.BASE_A == "GC") & \
            (nn_data.BASE_B == "init"),]
    
    if state < 20:  # if there is DNA duplex
        DNA_init = DNA_init + float(DNA_end_data.dH - temperature * DNA_end_data.dS) 
        if wt_target[state] == "A" or wt_target[state] == "T":
            DNA_start_data = nn_data.loc[(nn_data.NA_1 == "DNA") & \
                (nn_data.NA_2 == "DNA") & \
                (nn_data.BASE_A == "AT") & \
                (nn_data.BASE_B == "init"),]
        else:
            DNA_start_data = nn_data.loc[(nn_data.NA_1 == "DNA") & \
                (nn_data.NA_2 == "DNA") & \
                (nn_data.BASE_A == "GC") & \
                (nn_data.BASE_B == "init"),]
        DNA_init = DNA_init + float(DNA_start_data.dH - temperature * DNA_start_data.dS)
    return(RNA_base_stack + DNA_base_stack + DNA_init + RNA_init + PAM_dG)

def get_v_f_dG(state, dGs):
    if state==20:
        return(0)
    else:
        return(exp(-(dGs[state] - dGs[state - 1])/(2*R*temperature)))

def get_v_r_dG(state, dGs):
    if state==0:
        return(0)
    else:
        return(exp(-(dGs[state-1] - dGs[state])/(2*R*temperature)))

def get_mutant_dGs(wt_target, ref_dGs, Ts, Tv, mut_positions, base_lookup):
    mut_dGs = ref_dGs[:]
    for i in mut_positions:
        mut_dGs[i:] = [x if base_lookup[i] == wt_target[i - 1] else \
            x + Ts if ((wt_target[i - 1] in "AG" and base_lookup[i] in "AG") or (wt_target[i - 1] in "CT" and base_lookup[i] in "CT")) else \
            x + Tv if ((wt_target[i - 1] in "AG" and base_lookup[i] in "CT") or (wt_target[i - 1] in "CT" and base_lookup[i] in "AG")) else \
            np.nan for x in mut_dGs[i:]]
    return(mut_dGs)

def simulate_cluster_association(dna, max_time, halt_fraction, args):
    frac = 0
    cas_heap = []
    time = 0
    dna_molecules = [deepcopy(dna) for i in range(args.cluster_size)]
    for cluster_index in range(args.cluster_size):
        heappush(cas_heap, (dna_molecules[cluster_index].prepare_move(), cluster_index))

    while(time < max_time and frac < halt_fraction): # start simulation
        time, cluster_index = heappop(cas_heap)
        time_step = dna_molecules[cluster_index].move_and_prepare()
        frac = frac + dna_molecules[cluster_index].delta / args.cluster_size
        heappush(cas_heap, (time + time_step, cluster_index))

    return((time, dna_molecules))

def simulate_cluster_dissociation(on_dna_molecules, max_time, halt_fraction, args):
    off_dna_molecules = [deepcopy(dna) for dna in on_dna_molecules]
    for cluster_index in range(args.cluster_size):
        off_dna_molecules[cluster_index].v_fs[0] = 0

    frac = (args.cluster_size - sum([1 if dna.state > 0 else 0 for dna in off_dna_molecules])) / float(args.cluster_size)
    cas_heap = []
    time = 0
    for cluster_index in range(args.cluster_size):
        heappush(cas_heap, (off_dna_molecules[cluster_index].time_step, cluster_index))

    while(time < max_time and frac > 0. and frac > halt_fraction): # start simulation
        time, cluster_index = heappop(cas_heap)
        time_step = off_dna_molecules[cluster_index].move_and_prepare()
        frac = frac + off_dna_molecules[cluster_index].delta / args.cluster_size
        heappush(cas_heap, (time + time_step, cluster_index))

    return((time,off_dna_molecules))


# wt_time
# mut_state
# mut_base
# dna object

def simulate_sequence_association(param, args):
    single_dna = param[0]
    run_time = param[1]
    
    simulations = []
    halt_times = []

    for n in range(args.n_simulations):
        halt_time, dna_molecules = simulate_cluster_association(single_dna, run_time, 1., args)
        state_counts = [0 for i in range(21)]
        for cluster_index in range(args.cluster_size):
            state_counts[dna_molecules[cluster_index].state] = state_counts[dna_molecules[cluster_index].state] + 1

        simulations.append(state_counts)
        halt_times.append(halt_time)
    
    state_data = pd.DataFrame(simulations, columns = ["s" + str(i) for i in range(21)])
    mean_occupancy = (args.cluster_size - np.mean(state_data.s0)) / args.cluster_size
    mean_halt_time = np.mean(halt_times)
    return({'mutations' : ":".join([str(20 - s) + b for s, b in zip(single_dna.mut_states[::-1], single_dna.mut_bases[::-1])]), 'mut_states' : ",".join([str(s) for s in single_dna.mut_states]), 'mut_bases' : ",".join(single_dna.mut_bases), 'mean_halt_time' : mean_halt_time, 'mean_occupancy' : mean_occupancy})

def simulate_sequence_dissociation(param, args):
    single_dna = param[0]
    run_time = param[1]

    simulations = []
    halt_times = []
    
    for n in range(args.n_simulations):
        on_halt_time, on_dna_molecules = simulate_cluster_association(single_dna, run_time * args.time_multiplier, 1., args)
        off_halt_time, off_dna_molecules = simulate_cluster_dissociation(on_dna_molecules, run_time, 0., args)
        off_state_counts = [0 for i in range(21)]
        for cluster_index in range(args.cluster_size):
            off_state_counts[off_dna_molecules[cluster_index].state] = off_state_counts[off_dna_molecules[cluster_index].state] + 1

        simulations.append(off_state_counts)
        halt_times.append(off_halt_time)
    
    state_data = pd.DataFrame(simulations, columns = ["s" + str(i) for i in range(21)])
    mean_occupancy = (args.cluster_size - np.mean(state_data.s0)) / args.cluster_size
    mean_halt_time = np.mean(halt_times)
    return({'mutations' : ":".join([str(20 - s) + b for s, b in zip(single_dna.mut_states[::-1], single_dna.mut_bases[::-1])]), 'mut_states' : ",".join([str(s) for s in single_dna.mut_states]), 'mut_bases' : ",".join(single_dna.mut_bases), 'mean_halt_time' : mean_halt_time, 'mean_occupancy' : mean_occupancy})

def calc_single_mm_dna(wt_target, args):
    dna_info = []    
    dGs = [get_delta_G(i, [], {}, args.pam_binding_energy) for i in range(0,21)]
    dGs[args.adjustment_position:] = [x + args.energy_adjustment for x in dGs[args.adjustment_position:]]
    for mut_state in range(1, 21):
        for mut_base in "ACGT":
            if wt_target[mut_state - 1] == rc(mut_base):
                continue 
            dna_variant = dna(wt_target, dGs, [mut_state], [mut_base], args)
            dna_info.append(dna_variant)

    return(dna_info)

def calc_double_mm_dna(wt_target, args):
    dna_info = []
    dGs = [get_delta_G(i, [], {}, args.pam_binding_energy) for i in range(0,21)]
    dGs[args.adjustment_position:] = [x + args.energy_adjustment for x in dGs[args.adjustment_position:]]
    for first_state in range(1, 20):
        if first_state < args.second_state:
            for mut1_base in "ACGT":
                for mut2_base in "ACGT":
                    if wt_target[first_state - 1] == rc(mut1_base) or wt_target[args.second_state - 1] == rc(mut2_base):
                        continue 
                    dna_variant = dna(wt_target, dGs, [first_state, args.second_state], [mut1_base, mut2_base], args)
                    dna_info.append(dna_variant)

    return(dna_info)

if args.exit:
    sys.exit()



wt_dGs = [get_delta_G(i, [], {}, args.pam_binding_energy) for i in range(0,21)]
wt_dGs[args.adjustment_position:] = [x + args.energy_adjustment for x in wt_dGs[args.adjustment_position:]]
wt_dna = dna(wt_target, wt_dGs, [], [], args)

wt_halt_times = []
for i in range(args.n_wt_simulations):
    wt_time, wt_dnas = simulate_cluster_association(wt_dna, args.length_of_max_time, args.halt_fraction, args)
    wt_halt_times.append(wt_time)

mean_wt_time = np.mean(wt_halt_times)
dna_data = calc_single_mm_dna(wt_target, args)

if args.length_of_fixed_time:
    run_time = args.length_of_fixed_time
else:
    run_time = mean_wt_time

if args.dissociation:
    results = (Parallel(n_jobs=args.n_jobs)(delayed(simulate_sequence_dissociation)((x, run_time), args) for x in dna_data))
else:
    results = (Parallel(n_jobs=args.n_jobs)(delayed(simulate_sequence_association)((x, run_time), args) for x in dna_data))

run_data = pd.DataFrame(results)
run_data["mut_position_oldform"] = [20 - int(x) for x in run_data["mut_states"]]
run_data["mutations"] = [str(p) + b for p, b in zip(run_data["mut_position_oldform"], run_data["mut_bases"])]
merged_data = pd.merge(run_data, slope_data)
objective_value = euclidean(merged_data["slope"] / args.wt_slope , merged_data["mean_occupancy"] / args.halt_fraction)



summary_data = pd.DataFrame({\
    "filename" : (args.out_prefix).split("/")[-1] + ".all_sm.run.txt", \
    "energy_adjustment" :  args.energy_adjustment, \
    "adjustment_position" : args.adjustment_position, \
    "halt_fraction" : args.halt_fraction, \
    "transition_penalty" : args.transition_penalty, \
    "transversion_penalty" : args.transversion_penalty, \
    "association_rate" :  args.association_rate , \
    "euclidean" : objective_value}, \
    index=[0])


summary_data.to_csv(args.out_prefix + ".run_parameters.txt", sep="\t",index=False, na_rep="NA")
run_data.to_csv(args.out_prefix + ".all_sm.run.txt", sep="\t",index=False, na_rep="NA")

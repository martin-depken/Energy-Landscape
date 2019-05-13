import sys

PATH_HPC05 = '/home/svandersmagt/Energy_Landscape_dCas9/'
sys.path.append(PATH_HPC05)

import numpy as np
import SimulatedAnnealing_Nucleaseq_parallel as SimA
import calculate_cleavage_rate as clv
import functools
import Nucleaseq_data_processing as processing

## First set up a SA instance

filename = 'ECas9_cleavage_rate_and_y0_Canonical_OT-r_0-2.csv'
path_to_dataOn = PATH_HPC05 + 'data_nucleaseq_Finkelsteinlab/targetE/'
path_to_dataClv = PATH_HPC05 + 'data_nucleaseq_Finkelsteinlab/targetE/'

xdata, ydata, yerr = processing.prepare_multiprocessing_combined('1',filename,path_to_dataOn,path_to_dataClv,True)

perfectClv = np.float(len(ydata[0][0]))
perfectOn = np.float(len(ydata[0][1]))
singleClv = 0.0
singleOn = 0.0
doubleClv = 0.0
doubleOn = 0.0
for i in range(len(xdata)):
    if len(xdata[i])==1:
        singleClv += len(ydata[i][0])
        singleOn += len(ydata[i][1])
    if len(xdata[i])==2:
        doubleClv += len(ydata[i][0])
        doubleOn += len(ydata[i][1])
    
chi_weights = [1/perfectClv,1/singleClv,1/doubleClv,1/perfectOn,1/singleOn,1/doubleOn]

model = ['Clv_Saturated_general_energies_v2','general_energies_no_kPR']
Tstart = 1000
delta = 1.0
tol = 10^-5
Tfinal = 0.0
adjust_factor = 1.1
cooling_rate = 0.99
N_int = 1000
AR_low = 40.
AR_high = 60.
use_multiprocessing = True
nprocs = 19
use_relative_steps = False
objective_function = functools.partial(clv.calc_chi_squared,
                        guide_length=20,
                        model_id=model)
NMAC = False
reanneal = False
ReannealFactor = 1

upbnd = [10.0] + [10.0]*40 + [6.0] + [6.0] + [6.0]
lwrbnd = [0.0] + [-10.0]*20 + [0.0]*20 + [-1.0] + [1.0] + [3.0]
initial_guess =  [5.0] + [0.0]*20 + [5.0]*20 + [1.0] + [3.0] + [4.5]

SA = SimA.SimAnneal(model=model,
                   Tstart=Tstart,
                   delta=delta,
                   tol=tol,
                   Tfinal=Tfinal,
                    potential_threshold = np.inf,
                   adjust_factor=adjust_factor,
                   cooling_rate_high=cooling_rate,
                cooling_rate_low=cooling_rate,
                   N_int=N_int,
                    Ttransition=np.inf,
                   AR_low=AR_low,
                   AR_high=AR_high,
                   use_multiprocessing=use_multiprocessing,
                   nprocs=nprocs,
                   use_relative_steps=use_relative_steps,
                   objective_function=objective_function,
                   chi_weights=chi_weights,
                   NMAC=NMAC,
                   reanneal=reanneal)

it = 2000 #iterations per temperature
T = np.logspace(5,-1,200) #temperatures
Potentials = np.zeros([len(T),it/10]) #record potential every 10 iterations

output_file_results = PATH_HPC05 + 'TransitionT/result.txt'
Output = open(output_file_results,'w',1)  #third argument will force the I/O to write into the file every line
X = initial_guess

for j in range(len(T)):
    SA.potential = SimA.V(SA, xdata,ydata,yerr,X)
    for i in range(it):
        Xtrial = SimA.TakeStep(SA, X, lwrbnd, upbnd)
        Vnew = SimA.V(SA, xdata, ydata, yerr, Xtrial)
        Vold = SA.potential
        if (np.random.uniform() < np.exp(-(Vnew - Vold) / T[j])):
            X = Xtrial
            SA.potential = Vnew
        if i%10==0:
            Potentials[j][i/10] = SA.potential
            Output.write(str(SA.potential)+'\t')
    Output.write('\n')

Output.close()
for i in range(SA.nprocs):
    SA.inQ.put(None)
for w in SA.processes:
    w.join()
    
Var = np.var(Potentials,1)
C = Var/T**2
print C
'''
Code to simulate FRET experiments from energy landscape (d)Cas9
'''


import numpy as np
from scipy import linalg
import sys
sys.path.append('../code_Boyle/')
import CRISPR_dCas9_binding_curve_Boyle as dCas9
reload(dCas9);


def dwelltime_distribution(parameters, model_id, mismatch_positions, guide_length=20,
                           logbins=True, bins_per_decade=10, tMIN = 0.0001, tMAX=10**(6),
                           temporal_resolution=0.1):
    '''
        Generates dwelltime distribution from solution to Master Equation(s)
    :param parameters:
    :param model_id:
    :param mismatch_positions:
    :param guide_length:
    :param logbins:
    :param bins_per_decade:
    :param tMIN:
    :param tMAX:
    :param temporal_resolution:
    :return:
    '''

    # -- 1. Construct Master Equations using the parameter set (use function as before) ----
    rate_matrix = dCas9.get_master_equation(parameters, mismatch_positions, model_id, guide_length)

    # -- 2. Set concentration to zero: Single-molecule experiment ----
    new_rate_matrix = rate_matrix.copy()
    new_rate_matrix[0][0] = 0.0
    new_rate_matrix[1][0] = 0.0

    # -- 3. initial condition: Molecule enteres at PAM state (assume timeres is good enough) ----
    PAM_bound = np.zeros(guide_length + 2)
    PAM_bound[1] = 1.0

    # -- 4. construct the bins/ timepoints at which we evaluate the distribution ----


    # Log-spaced bins:
    if  logbins:
        times = [tMIN]
        t = tMIN
        while t <= tMAX:
            t = times[-1] * 10 ** ((1 / float(bins_per_decade)))
            times.append(t)
    else:
        # Linear bins:
        dt = temporal_resolution
        times = np.arange(tMIN, tMAX, dt)

    # -- 5. For different timepoints, solve for the dwelltime distribution ----
    dwelltime_dist = []
    Probabilities = PAM_bound
    for i in range(len(times) - 1):
        dt = times[i + 1] - times[i]
        matrix_exponent = linalg.expm(+new_rate_matrix * dt)
        Probabilities = matrix_exponent.dot(Probabilities)
        rate_to_sol = np.diag(new_rate_matrix, k=1)[0]
        dwelltime_dist.append(Probabilities[1] * rate_to_sol)
    return dwelltime_dist, times[1:]


def simulate_FRET_trace(parameters, model_id,
                        mismatch_positions,
                        filename,
                        concentration=10.,
                        guide_length=20,
                        measurement_time=300.0,
                        temporal_resolution=0.1,
                        noise_amplitude=0.0):
    '''
    Use Gillespie simulation to generate a mock FRET trace that takes into account the Cas9 (or DNA) concentration,
    finite time resulotion of camera, finite measurement time (possibly due to photobleaching) and experimental
    noise in the FRET values (resulting from noise in the fluoresence levels)

    :param parameters: from SA fit
    :param model_id:
    :param mismatch_positions:
    :param filename:
    :param concentration:
    :param guide_length:
    :param measurement_time:
    :param temporal_resolution:
    :param noise_amplitude:
    :return: writes the trace to a file
    '''
    # ------- use the parameters from fit to determine rates in free-energy landscape -----
    epsilon, forward_rates = dCas9.unpack_parameters(parameters, model_id, guide_length)

    # ------ adjust to concentration of experiment: Boyle parameters are at 10nM ------
    epsilon[0] -= np.log(concentration / 10.0)
    forward_rates[0] *= concentration / 10.0

    # ------ deduce the free-energy landscape for the construct at hand -------
    energies = dCas9.get_energies(epsilon, mismatch_positions, guide_length)
    backward_rates = dCas9.get_backward_rates(energies, forward_rates, guide_length)

    # ------ perform gillespie algorithm to mock a FRET signal -------------------
    ###  states:
    # -1: unbound
    #  0: PAM,
    #  1-20: number of bp in R-loop

    # At t=0 start in solution:
    state = -1
    time = 0.0

    # Set (average) or ideal FRET values for bound and unbound configurations:
    high_FRET = 0.8
    low_FRET = 0.0  # maybe it is actually not zero?

    O = open(filename, 'w')
    while time < measurement_time:

        # ----- Gillespie Algorithm: Determine when next switch in FRET -----
        rate_fwd = forward_rates[state + 1]
        rate_back = backward_rates[state + 1]
        rate = rate_fwd + rate_back
        dwell_time = np.random.exponential(scale=rate ** (-1))

        # ---- Adjust the FRET value ----
        if state == -1:
            FRET = low_FRET
        else:
            FRET = high_FRET

        # ----  Current event starts at:
        time_last_event = time

        # ---- make a time trace with finite sampling resolution ----
        while time <= (time_last_event + dwell_time):
            time += temporal_resolution

            # ---- Add noise to FRET value ----
            FRET += np.random.normal(loc=0.0, scale=noise_amplitude)

            # ---- FRET values cannot fall below zero (negative fluoresence) ---
            FRET = np.maximum(FRET, 0.0)

            # ---- store the datapoint ----
            # ---- hold the (idealised) FRET value between two consequetive events ----
            O.write(str(time) + '\t')
            O.write(str(FRET) + '\t')
            O.write(str(state) + '\n')

        # ---- Now use Gillespie Algorithm to determine the nature of the next event
        U = np.random.uniform()
        if U < rate_fwd / (rate_fwd + rate_back):
            state += 1
        else:
            state -= 1

    O.close()
    return
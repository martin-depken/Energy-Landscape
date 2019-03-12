import numpy as np

def get_backward_rates(energies, forwardrates, guide_length=20):
	# kb(n) = kf(n-1) * exp(+epsilon(n))
    backwardrates = np.zeros(guide_length)
    backwardrates = forwardrates[:-1] * np.exp(energies)
    #backwardrates = forwardrates[:-1] * np.exp(energies) # version from lab
    #backwardrates[:-1] = forwardrates[:-1] * np.exp(energies[:-1]) # XXX CHECK THIS!
    return backwardrates
	
def build_rate_matrix(forward_rates, backward_rates):
    diagonal1 = - (forward_rates + backward_rates)
    diagonal2 = backward_rates[1:]
    diagonal3 = forward_rates[:-1]
    rate_matrix = np.diag(diagonal1, k=0) + np.diag(diagonal2, k=1) + np.diag(diagonal3, k=-1)
    return rate_matrix

def mean_first_passage_time(M):
	guidelength = len(M)
	onevec = np.ones(len(M))
	everything_unbound = np.array([1.0] + [0.0] * (guidelength-1))
	return onevec.dot(np.linalg.solve(-M,everything_unbound))
	
def get_backward_rates_2(energies, forward_rates,guide_length=20):
    '''
    Apply detailed balance condition to get the backward rates from the energies and forward rates

    :param energies:
    :param forward_rates:
    :param guide_length:
    :return:
    '''
    # 0) Construct array containing backward rates
    backward_rates = np.zeros(guide_length+1)

    # 1) Apply detailed balance condition:
    backward_rates = forward_rates[:-1] * np.exp(energies)
    return backward_rates

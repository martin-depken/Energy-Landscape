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
	onevec = np.ones(len(M))
	everything_unbound = np.array([1.0] + [0.0] * (guidelength-1))
	return onevec.dot(np.linalg.solve(-M,everything_unbound))
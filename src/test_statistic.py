from numpy import dot, cov, array, mean, std, arange, sqrt, subtract
from numpy.linalg import inv

def get_amplitudes(cloud, observation, modes):
	cloud_amps = dot(cloud, modes)
	observed_amps = dot(observation, modes)

	return cloud_amps, observed_amps


def distance_modulus(cloud_amps, observed_amps):

	covariance = cov(cloud_amps.T)
	means = mean(cloud_amps, axis = 0)

	try:
		distance_modulus = dot((observed_amps - means).T , dot(inv(covariance), (observed_amps - means)))
	except:
		distance_modulus = (observed_amps - means)**2/covariance
	return distance_modulus


def global_statistic(cloud_amps, observed_amps):
	numsamples, numparams = cloud_amps.shape
	cloud_distances = array([distance_modulus(cloud_amps, cloud_amps[i,:]) for i in arange(numsamples)])
	observed_distance = distance_modulus(cloud_amps, observed_amps)
	expectation = mean(cloud_distances)
	sigma = std(cloud_distances)
	chi2 = (observed_distance-expectation)**2/sigma**2
	expected_chi2 = subtract(cloud_distances, expectation)**2/sigma**2
	zscore = (chi2-mean(expected_chi2))/std(expected_chi2)
	return chi2, zscore


def max_statistic(cloud_amps, observed_amps):

	chi2 = -1.0e-30
	zscore = -1.0e-30
	numsamples, numparams = cloud_amps.shape
	for i in arange(1,numparams):
		cloud_distances = array([distance_modulus(cloud_amps[:,:i], cloud_amps[j, :i]) for j in arange(numsamples)])
		observed_distance = distance_modulus(cloud_amps[:,:i], observed_amps[:i])

		expectation = mean(cloud_distances)
		sigma = std(cloud_distances)
		expected_chi2 = subtract(cloud_distances , expectation)**2/sigma**2
		newscore = (chi2-mean(expected_chi2))/std(expected_chi2)
		if abs(newscore) > abs(zscore):
			chi2 = (observed_distance - expectation)**2/sigma**2
			zscore = newscore
	return chi2, zscore





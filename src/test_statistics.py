import model
import data
import numpy as np

class test_statistic(dict):
	def __init__(self, compression_scheme, experiments, *args, **kwargs):
		super(compression, self).__init__(*args, **kwargs)
		self.scheme = compression_scheme
		self.data = self.experiments
		try:
			self.experiments[0]
		except:
			self.experiments = [self.experiments]

		self.nummodes = self.scheme.nummodes

	def get_measured_amplitudes(self):

		self.measured_amps = compute.get_measured_amps(self)

	def get_lcdm_cloud(self):

		self.lcdm_amps = compute.get_lcdm_amps(self)
		self.expected_amps = np.mean(self.lcdm_amps, axis = 0)
		self.amp_covariance = np.cov(self.lcdm_amps)

	def global_stat(self):

		self.get_measured_amplitudes()
		self.get_lcdm_cloud()
		dx = self.measured_amps - self.expected_amps
		self.global_chi2 = np.dot(dx, np.dot(self.amp_covariance, dx))

		self.global_zscore = dx / diag(self.amp_covariance)

	def max_stat(self):

		

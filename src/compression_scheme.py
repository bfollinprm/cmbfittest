import model
import data
import numpy as np
import compute


class compression_scheme(dict):
	def __init__(self, compression_type, *args, **kwargs):
		super(compression_scheme, self).__init__(*args, **kwargs)
		self.__dict__ = self

		self.type = compression_type

		try:
			self.dcldphi
		except:
			self.dcldphi = np.identity(5001)

		try:
			self.noise_cov
		except:
			self.noise_cov = np.identity(5001)
		try:
			self.supermodel_cov
		except:
			self.supermodel_cov = np.outer(self.dcldphi, self.dcldphi)

		try:
			self.SM_cov
		except:
			self.SM_cov = np.identity(5001)

		try:
			self.kappa
		except:
			self.kappa = 0

		try:
			self.nummodes
		except:
			self.nummodes = 30


	def get_compression_modes(self):
		self.eigspectrum, self.modes = compute.get_modes(self)

		self.modes = self.modes[:nummodes]

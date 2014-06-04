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
			self.dcldphi = np.identity(2501)

		try:
			self.noise_cov
		except:
			self.noise_cov = np.identity(2501)
		try:
			self.supermodel_cov
		except:
			self.supermodel_cov = np.dot(self.dcldphi.T, self.dcldphi)

		try:
			self.SM_cov
		except:
			self.SM_cov = np.identity(2501)

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

		self.modes = self.modes[:,:self.nummodes]


	def save_modes(self):
		root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
		filename = root + '/data/compression_modes/%s_modes.npy'%self.type

		try:
			numpy.save(filename, self.modes)
		except:
			self.dcldphi = self.get_compression_modes()
			numpy.save(filename, self.modes)

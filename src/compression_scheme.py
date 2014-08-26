import model
import data
import numpy as np
from numpy.linalg import eig, inv
from scipy.linalg import block_diag
from pylab import *


class compression_scheme(dict):
	def __init__(self, compression_type, *args, **kwargs):
		super(compression_scheme, self).__init__(*args, **kwargs)
		self.__dict__ = self

		self.type = compression_type

		try:
			self.ellmax
		except:
			self.ellmax = 2500

		try:
			self.dcldphi
		except:
			self.dcldphi = np.identity(self.ellmax + 1)

		try:
			self.noise_cov
		except:
			self.noise_cov = [np.identity(self.ellmax + 1)]
		try:
			self.supermodel_cov
		except:
			try:
				self.supermodel_cov = np.dot(self.dcldphi.T, self.dcldphi)
			except:
				self.supermodel_cov = np.identity(self.ellmax + 1)
		try:
			self.SM_cov
		except:
			self.SM_cov = np.identity(self.ellmax + 1)

		try:
			self.kappa
		except:
			self.kappa = 0

		try:
			self.nummodes
		except:
			self.nummodes = 30
		try:
			self.windows
		except:
			self.windows = [np.identity(self.ellmax + 1)]

	def get_compression_modes(self):
		self.eigspectrum, self.modes = get_modes(self)

		self.modes = self.modes[:,:self.nummodes]

	def save_modes(self):
		root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
		filename = root + '/data/compression_modes/%s_modes.npy'%self.type

		try:
			numpy.save(filename, self.modes)
		except:
			self.dcldphi = self.get_compression_modes()
			numpy.save(filename, self.modes)


##### Auxiliary Functions #####
def get_modes(test):
	if test.type in ['LDA', 'PCA', 'S_NL', 'S_N']:
		num_experiments = len(test.noise_cov)
		SM_cov = []
		supermodel_cov = []
		off_diag = []
		noise_cov = []
		for i in np.arange(num_experiments):
			SM_cov.append(
							np.dot(
								np.dot(test.windows[i], test.SM_cov),
								test.windows[i].T
								  )
						 )
			supermodel_cov.append(
							np.dot(
								np.dot(test.windows[i], test.supermodel_cov),
								test.windows[i].T
								  )
						 )
			noise_cov.append(test.noise_cov[i])

		supermodel_cov = block_diag(*supermodel_cov)
		SM_cov = block_diag(*SM_cov)
		### Add the off diagonal terms
		for i in np.arange(num_experiments):
			for j in np.arange(i):
				xstart = sum([test.windows[k].shape[0] for k in arange(i)])
				ystart = sum([test.windows[k].shape[0] for k in arange(j)])
				SM_cov[xstart:test.windows[i].shape[0] + xstart, ystart:test.windows[j].shape[0] + ystart] = dot(test.windows[i], dot(test.SM_cov, test.windows[j].T))
				supermodel_cov[xstart:test.windows[i].shape[0] + xstart, ystart:test.windows[j].shape[0] + ystart] = dot(test.windows[i], dot(test.supermodel_cov, test.windows[j].T))
				SM_cov[ystart:test.windows[j].shape[0] + ystart, xstart:test.windows[i].shape[0] + xstart] = dot(test.windows[i], dot(test.SM_cov, test.windows[j].T)).T
				supermodel_cov[ystart:test.windows[j].shape[0] + ystart, xstart:test.windows[i].shape[0] + xstart] = dot(test.windows[i], dot(test.supermodel_cov,test.windows[j].T)).T

		figure()
		matshow(log(abs(SM_cov[48:,48:])))
		colorbar()
		figure()
		matshow(log(abs(supermodel_cov)))
		colorbar()
		show()
		noise_cov = block_diag(*noise_cov)
		print test.kappa
		# 	#print experiment_test_matrices[-1].shape
		w, v = eigh(noise_cov + test.kappa * SM_cov)
		denominator = dot(v, dot(diag(1/w), v.T))
		test_matrix = np.dot(denominator, supermodel_cov)
		# figure()
		# matshow(log(abs(test_matrix[46:,46: ])))
		# colorbar()
		# show()
		# figure()
		matshow(log(abs(test_matrix)))
		colorbar()
		show()
		# figure()
		# loglog(abs(real(eigvals(test_matrix[:45,:45])[::-1])))
		# loglog(abs(imag(eigvals(test_matrix[:45,:45])[::-1])))
		# show()
		# figure()
		# loglog(abs(real(eigvals(test_matrix[46:,46:])[::-1])))
		# loglog(abs(imag(eigvals(test_matrix[46:,46:])[::-1])))

		# show()
		w, modes = eig(test_matrix)

	
	# if test.type in ['fisher', 'FN', 'FNL']:
	# 	test_matrix = np.dot(inv(test.noise_cov[2:,2:]) + test.kappa + test.SM_cov[2:,2:])
	# 	w, modes = eig(test_matrix)

	# padded_modes = np.zeros([test.ellmax + 1, test.nummodes])
	# for i in np.arange(w.size):
	# 	for j in np.arange(num_experiments):
	# 		location = sum([test.noise_cov[k].shape[0] for k in np.arange(j)]+2 * j)
	# 		np.insert(np.real(modes[:,i]), location, np.zeros(2))

	#Make sure they're sorted
	ii = np.argsort(w)
	w = w[ii]
	w = w[::-1]
	modes = modes[:,ii]
	modes = modes[:,::-1]
	figure()
	loglog(w)
	show()
	return np.real(w), modes 

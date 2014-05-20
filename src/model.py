import numpy
import os.path as osp
import compute

class model(dict):
	'''
	class to store a model, which is (a) a set of linearized responses to parameters dcldphi, and (b) a covariance matrix
	'''
	def __init__(self, modelspace, *args, **kwargs):
		super(model, self).__init__(*args, **kwargs)
		self.__dict__ = self
		self.model = modelspace

		try:
			self.force_calculation
		except:
			self.force_calculation = False


	def get_dcldphi(self):
		'''
		Computes the linearized responses in C_l space from varying the parameters of the model, or grabs them from a 
		previous run if they're saved.
		'''
		root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
		filename = root + '/data/cl_fid.npz'

		try:
			fiducial_model = numpy.load(filename)
			self.fiducial_cl = fiducial_model['cl_TT']
			self.fiducial_params = fiducial_model['params'].item()
		except:
			raise Exception('Critical Error: No fiducial model found, or fiducial model in incorrect format.  Need file %s to contain dictionary with keys cl_TT and params'%filename)


		if self.force_calculation == False:
			## Try to load the data, if not revert to the calculation
			root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
			filename = root + '/data/models/%s_dcldphi.npy'%self.model
			try:
				dcldphi = numpy.load(filename)
			except:
				print 'Error: file %s not found.  Calculating dcldphi'%filename
				dcldphi = compute.calc_dcldphi(self)
		else:
			## Do the calculation
			dcldphi = compute.calc_dcldphi(self)

		self.dcldphi  = dcldphi

	def get_covariance(self):
		'''
		Calculates the covariance matrix of the model under some prior distribution of the parameters.
		'''


		### Make sure the modes are loaded
		try:
			self.dcldphi
		except:
			self.get_dcldphi()

		try:
			self.prior
		except:
			self.get_default_prior()

		###And project back into multipole space
		self.covariance = numpy.dot(self.dcldphi, numpy.dot(self.prior, self.dcldphi.T))

	def get_default_prior(self):
		'''
		provides a default covariance structure on the modes; currently the identity matrix
		'''
		self.prior = numpy.identity(self.dcldphi.shape[1])

	def save_dcldphi(self):
		root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
		filename = root + '/data/models/%s_dcldphi.npy'%self.model

		try:
			numpy.save(filename, self.dcldphi)
		except:
			self.dcldphi = self.get_dcldphi()
			numpy.save(filename, self.dcldphi)









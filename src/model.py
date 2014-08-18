
import numpy as np
import os.path as osp
import os.path as osp
from numpy import pi
from numpy.random import multivariate_normal as normal
from numpy.linalg import eig, inv

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


		try:
			self.ellmax
		except:
			self.ellmax = 5000


	def get_dcldphi(self):
		'''
		Computes the linearized responses in C_l space from varying the parameters of the model, or grabs them from a 
		previous run if they're saved.
		'''
		root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
		filename = root + '/data/cl_fid.npz'

		try:
			fiducial_model = np.load(filename)
			self.fiducial_cl = fiducial_model['cl_TT'][:self.ellmax + 1]
			self.fiducial_params = fiducial_model['params'].item()
		except:
			raise Exception('Critical Error: No fiducial model found, or fiducial model in incorrect format.  Need file %s to contain dictionary with keys cl_TT and params'%filename)


		if self.force_calculation == False:
			## Try to load the data, if not revert to the calculation
			root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
			filename = root + '/data/models/%s_dcldphi.npy'%self.model
			try:
				dcldphi = np.load(filename)
			except:
				print 'Error: file %s not found.  Calculating dcldphi'%filename
				dcldphi = calc_dcldphi(self)
		else:
			## Do the calculation
			dcldphi = calc_dcldphi(self)

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
		self.covariance = np.dot(self.dcldphi, np.dot(self.prior, self.dcldphi.T))

	def get_default_prior(self):
		'''
		provides a default covariance structure on the modes; currently the identity matrix
		'''
		self.prior = np.identity(self.dcldphi.shape[1])

	def save_dcldphi(self):
		root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
		filename = root + '/data/models/%s_dcldphi.npy'%self.model

		try:
			np.save(filename, self.dcldphi)
		except:
			self.dcldphi = self.get_dcldphi()
			np.save(filename, self.dcldphi)


##### Auxiliary Functions #####

def calc_dcldphi(model):
	
	## Get the Boltzmann code algorithm to calculate responses in C_l from changing parameters phi
	try:
		import classy
	except:
		raise Exception('Critical Error: classy from www.class-code.net needed to calculate power spectra.  Not found')



	### Load the fiducial model and fiducial parameters, if they exist.
	fiducial_params = model.fiducial_params
	fiducial_cl = model.fiducial_cl[:model.ellmax + 1]
	observables = classy.Class()
	ell = np.arange(model.ellmax + 1)

	if 'lcdm' == model.model:
		parameter_shifts = {'A_s':5.0e-10, 'n_s':0.005, 'tau_reio':0.005, 'omega_b':0.00001, 'omega_cdm':0.0001, 'H0':2}
		dcldphi = np.zeros([model.ellmax + 1, len(parameter_shifts)])
		for i, param in enumerate(parameter_shifts.keys()):
			observables.set(**fiducial_params)
			observables.set(**{param:parameter_shifts[param] + fiducial_params[param]})
			print 'Finding response for varying %s'%param
			try:
				observables.compute()
			except:
				raise Exception('Critical Error: unable to compute spectra--likely bad parameters in fiducial model.  See class error.')
			try: 
				cl = observables.lensed_cl(lmax = max([5000, model.ellmax + 100]))['tt'][:model.ellmax + 1]
			except:
				raise Exception('Critical Error:  power spectrum not defined out to ell = 5000--likely bad parameter in fiducial model.  See class error.')

			dcldphi[:,i] =  (cl * ell * (ell + 1)/(2 * pi) * (2.73 * 1.0e6)**2 - fiducial_cl)/parameter_shifts[param]
			observables.cleanup()

	if 'primordial' == model.model:
		perturbation_locations = np.arange(-8, 2, 0.05)
		dcldphi = np.zeros([model.ellmax + 1, perturbation_locations.size])
		for i, lnk in enumerate(perturbation_locations):
			observables.set(**fiducial_params)
			observables.set(perturb_loc = lnk)
			observables.set(perturb_height = 5.0e-11)
			observables.set(perturb_width = 0.05)
			print 'finding response for primordial perturbation centered at lnk = %g'%lnk
			try:
				observables.compute()
			except:
				observables.compute()
				raise Exception('Critical Error: unable to compute spectra--likely bad parameters in fiducial model.  See class error.')

			try: 
				cl = observables.lensed_cl(lmax = 5000)['tt'][:model.ellmax + 1]
			except:
				raise Exception('Critical Error:  power spectrum not defined out to ell = 5000--likely bad parameter in fiducial model.  See class error.')

			dcldphi[:,i] =  (cl * ell * (ell + 1)/(2 * pi) * (2.73 * 1.0e6)**2 - fiducial_cl)#/5.0e-11
			observables.cleanup()

	if 'transfer' == model.model:
		parameter_shifts = {
			'Omega_k':0.005, 
			'Omega_fld':1 - 1.0e4 * (fiducial_params['omega_b'] + fiducial_params['omega_cdm'])/fiducial_params['H0']**2, 
			'YHe':.26, 
			'N_ur':3.146,
			'N_ncdm':1
			}
		dcldphi = np.zeros([model.ellmax + 1, len(parameter_shifts)])
		for i, param in enumerate(parameter_shifts.keys()):
			observables.set(**fiducial_params)
			dphi = parameter_shifts[param]
			if param == 'N_ncdm':
				observables.set(N_ur = 2.046)
				observables.set(m_ncdm = 0.06)
				dphi = 0.06
			if param == 'Omega_fld':
				observables.set({'w0_fld':-0.95, 'wa_fld':0, 'cs2_fld':1})
				dphi = 0.05
			if param == 'N_ur':
				dphi = 0.1
			if param == 'YHe':
				dphi = 0.1
			observables.set(**{param:parameter_shifts[param]})

			try:
				observables.compute()
			except:
				raise Exception('Critical Error: unable to compute spectra--likely bad parameters in fiducial model.  See class error.')
			try: 
				cl = observables.lensed_cl(lmax = 5000)['tt']
			except:
				raise Exception('Critical Error:  power spectrum not defined out to ell = 5000--likely bad parameter in fiducial model.  See class error.')

			dcldphi[:,i] =  (cl * ell * (ell + 1)/(2 * pi) * (2.73 * 1.0e6)**2 - fiducial_cl)/dphi
			observables.cleanup()


	########### Add more model space definitions here. ############
	###															###
	###															###
	###															###
	###															###
	############################################################### 
	return dcldphi




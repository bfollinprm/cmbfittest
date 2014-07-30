import numpy as np
import os.path as osp
from numpy import pi
from numpy.random import multivariate_normal as normal
from numpy.linalg import eig, inv
from pylab import *

#import simulate


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
				print 
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


def calc_covariance(data):


	root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
	fileroot = root + '/data/experiments/%s/'%data.experiment

	if data.experiment == 'spt':
		filename = fileroot + 'covariance.dat'
		covariance = np.loadtxt(filename)
		

	if data.experiment == 'wmap':
		fidCl = data.fiducial_cl[:data.ellmax + 1]
		ondiag = np.loadtxt(fileroot + '/data/highl/wmap_likelihood_inputs_tt.p4v6.wmap9.kq85.cinv_v3.dat')
		offdiag = np.loadtxt(fileroot + '/data/highl/wmap_likelihood_inputs_tt_offdiag.p4v4.wmap9.kq85.cinv_v3.dat')


		ell = ondiag[:,0].astype('int')
		bandpower = ondiag[:,1]
		Nl = ondiag[:,2]
		Fsky = ondiag[:,3]
		l1 = offdiag[:,0].reshape([1199, 599]).astype('int')
		l2 = offdiag[:,1].reshape([1199, 599]).astype('int')
		Epsilon_temp = offdiag[:,2].reshape([1199, 599])
		R_eff_temp = offdiag[:,3].reshape([1199, 599])

		#R_eff = np.diag(np.zeros(data.ellmax + 1))
		R_eff = np.diag(np.zeros(1201))
		R_eff[l1, l2] = R_eff_temp
		R_eff[l2, l1] = R_eff_temp


		Epsilon = np.diag(np.zeros(1201))
		#np.diag(np.zeros(data.ellmax + 1))
		Epsilon[l1, l2] = Epsilon_temp
		Epsilon[l2, l1] = Epsilon_temp
		tttt = np.zeros(1201) + 1e-30
		#np.zeros(data.ellmax + 1) + 1
	
		tttt[ell] = 2 * (fidCl[ell] + Nl)**2/(2*ell + 1)/Fsky**2


		tttt = np.outer(np.sqrt(tttt), np.sqrt(tttt))

		inv_covariance = np.diag(1.0/np.diag(tttt)) - R_eff/np.sqrt(tttt) + Epsilon/tttt
		imshow(log(inv_covariance[2:,2:]))
		colorbar()
		show()
		covariance = np.diag(np.zeros(data.ellmax + 1)+1.0e30)
		covariance[2:1201, 2:1201] = inv(inv_covariance[2:,2:])*1e12


	if data.experiment == 'planck_like':
		ell = np.arange(data.ellmax + 1)
		Fsky = 0.3
		Nl  = 0.00015048 * np.exp(ell**2 * 7.69112e-7 ) * (ell**2) / 2.0 / np.pi
		covariance = 2 * (data.fiducial_cl[:data.ellmax + 1] + Nl)**2/(2*ell + 1)/Fsky
		covariance = np.diag(covariance)
	return covariance



def load_window_functions(data):

	root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
	fileroot = root + '/data/experiments/%s/'%data.experiment

	if data.experiment == 'spt':

		windows = [np.loadtxt(fileroot + 'windows/window_%i'%i) for i in range(1,48)]
		ell = np.array(windows)[0,:,0].astype('int')
		weights = np.array(windows)[:,:,1]
		windows = np.zeros([47,data.ellmax + 1])
		ell = [i for i in ell if i < data.ellmax + 1]
		windows[:, ell] = weights[:,np.array(ell)-50]
		return windows

	if data.experiment == 'wmap':

		windows = np.identity(data.ellmax + 1)[:1200,:]
		return windows

	if data.experiment == 'planck_like':

		windows = np.identity(data.ellmax + 1)
		return windows

	return windows


def get_bandpowers(data):
	root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
	fileroot = root + '/data/experiments/%s/'%data.experiment

	if data.experiment == 'spt':

		bandpowers = np.loadtxt(fileroot + 'bandpowers.dat')

	if data.experiment == 'wmap':

		ondiag = fileroot + 'wmap_likelihood_inputs_tt.p4v6.wmap9.kq85.cinv_v3.dat'
		ell = ondiag[:,0]
		bandpowers_tmp = ondiag[:,1]
		bandpowers = np.zeros(data.ellmax + 1)
		bandpowers[ell] = bandpowers_tmp


	if data.experiment == 'planck_like':

		bandpowers = data.fiducial_cl[:data.ellmax + 1] + normal(np.zeros(data.covariance.shape[0]), data.covariance)

	return bandpowers


def get_modes(test):

	if test.type in ['LDA', 'PCA', 'S_NL', 'S_N']:
		test_matrix = np.dot(inv(test.noise_cov[2:,2:] + test.kappa * test.SM_cov[2:,2:]), test.supermodel_cov[2:,2:])
		w, modes = eig(test_matrix)
	
	if test.type in ['fisher', 'FN', 'FNL']:
		test_matrix = np.dot(inv(test.noise_cov[2:,2:]) + test.kappa + test.SM_cov[2:,2:])
		w, modes = eig(test_matrix)

	padded_modes = np.zeros([test.ellmax + 1, test.ellmax - 1])
	for i in np.arange((np.real(modes).shape)[0]):
		padded_modes[2:,i] = np.real(modes[:,i])

	##Make sure they're sorted
	#ii = np.argsort(w)
	#w = w[ii]
	#w = w[::-1]
	#padded_modes = modes[:,ii]
	#padded_modes = modes[:,::-1]
	return np.real(w), np.real(padded_modes) 

	




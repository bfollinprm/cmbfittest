import numpy as np
import os.path as osp
import common
import os.path as osp
from numpy import pi
from numpy.random import multivariate_normal as normal
from numpy.linalg import eig, inv
from pylab import *

class data(dict):
	'''
	class to store an experiment, which is defined as a noise model, window functions, and a set of bandpowers.
	'''

	def __init__(self, experiment, *args, **kwargs):
		super(data, self).__init__(*args, **kwargs)
		self.__dict__ = self
		self.experiment = experiment

		try: 
			self.ellmax
		except:
			self.ellmax = 2500


	def get_covariance(self):
		self.lslice = get_lslice(self)
		root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
		filename = root + '/data/cl_fid.npz'
		try:
			fiducial_model = np.load(filename)
			self.fiducial_cl = fiducial_model['cl_TT']
		except:
			raise Exception('Critical Error: No fiducial model found, or fiducial model in incorrect format.  Need file %s to contain dictionary with keys cl_TT and params'%filename)

		self.covariance = calc_covariance(self)


	def get_window_functions(self):
		self.windows = load_window_functions(self)


	def get_bandpowers(self):
		try:
			self.observation
		except:

			self.observation = get_bandpowers(self)



###### Auxiliary functions ######


def get_lslice(data):
	if data.experiment == 'spt':
		lmax = min(data.ellmax+1, 3300)
		lslice = np.arange(50,lmax)
	if data.experiment == 'wmap':
		lmax = min(data.ellmax+1, 1201)

		lslice = np.arange(2,lmax)
	if data.experiment == 'planck_like':
		lmax = min(data.ellmax+1, 2501)

		lslice = np.arange(2,lmax)
	return lslice

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
		beam_info = np.loadtxt('./data/experiments/wmap/data/highl/top_ten_modes.beam_covariance_VW_combined.dat')

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
		tttt = np.zeros(1201)
		#np.zeros(data.ellmax + 1) + 1

		beam_modes = np.zeros([10,1201])
		beam_modes[:,ell] = beam_info[:,2].reshape(10,1199)
		for i in np.arange(9):
			beam_modes[i,ell] = fidCl[ell] * beam_modes[i, ell]
		beam_modes[9, ell] = .1 * beam_modes[9,ell]
	
		tttt[ell] = 2 * (fidCl[ell] + Nl)**2/(2*ell + 1)/Fsky**2

		tttt = np.outer(np.sqrt(tttt), np.sqrt(tttt))

		inv_covariance = np.diag(1.0/np.diag(tttt)) - R_eff/(tttt) + Epsilon/tttt**2
		covariance = np.diag(np.zeros(data.ellmax + 1))


		covariance[2:1201, 2:1201] = inv(inv_covariance[2:1201,2:1201])
		covariance = covariance[2:1201, 2:1201]
		# for i in np.arange(10):
		# 	covariance[2:1201, 2:1201] += np.outer(beam_modes[i,ell], beam_modes[i,ell])


	if data.experiment == 'planck_like':
		ell = np.arange(2,data.ellmax + 1)
		Fsky = 0.3
		Nl  = 0.00015048 * np.exp(ell**2 * 7.69112e-7 ) * (ell**2) / 2.0 / np.pi
		covariance = 2 * (data.fiducial_cl[2:data.ellmax + 1] + Nl)**2/(2*ell + 1)/Fsky
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
		windows[:, np.array(ell)] = weights[:,np.array(ell)-50]
		return windows

	if data.experiment == 'wmap':
		windows = np.identity(data.ellmax + 1)[2:1201,:]
		return windows

	if data.experiment == 'planck_like':

		windows = np.identity(data.ellmax + 1)[2:2501, :]
		return windows

	return windows


def get_bandpowers(data):
	root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
	fileroot = root + '/data/experiments/%s/'%data.experiment

	if data.experiment == 'spt':
		bandpowers = np.loadtxt(fileroot + 'bandpowers.dat')

	if data.experiment == 'wmap':

		ondiag = np.loadtxt(fileroot + 'data/highl/wmap_likelihood_inputs_tt.p4v6.wmap9.kq85.cinv_v3.dat')

		ell = ondiag[:,0].astype('int')
		bandpowers = ondiag[:,1]


	if data.experiment == 'planck_like':
		bandpowers = data.fiducial_cl[2:data.ellmax + 1] + normal(np.zeros(data.covariance.shape[0]), data.covariance)

	return bandpowers







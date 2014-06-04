import numpy as np
import os.path as osp

from numpy.random import multivariate_normal as normal
from numpy.linalg import inv
from cosmoslik import param_shortcut, get_plugin, SlikDict, SlikPlugin, Slik

from scipy.optimize import minimize



# param = param_shortcut('start','scale')


# class MCMC_chain(SlikPlugin):
# 	def __init__(self, experiments):
#         super(SlikPlugin,self).__init__()
#         self.__dict__ = self

#         self.experiments = experiments

# 		try: 
# 			self.experiments[0]
# 		except:
# 			self.experiments = [self.experiments]


#         self.cosmo = get_plugin('models.cosmology')(
#             logA = param(3.2),
#             ns = param(0.96),
#             ombh2 = param(0.0221),
#             omch2 = param(0.12),
#             tau = param(0.09,min=0,gaussian_prior=(0.0925,0.015)),
#             theta = param(0.010413),
#             omnuh2 = 0.000645,
#         )
                      


# 		fileroot, extension =  osp.split(osp.dirname(osp.abspath(__file__)))


#         self.get_cmb = get_plugin('models.pico')(
#             datafile=fileroot + '/data/pico3_tailmonty_v33.dat'
#         )

#         self.bbn = get_plugin('models.bbn_consistency')()
#         self.hubble_theta = get_plugin('models.hubble_theta')()  
#         self.priors = get_plugin('likelihoods.priors')(self)


#         self.sampler = get_plugin('samplers.metropolis_hastings')(
#             self,
# 	     	num_samples=10000,
#             proposal_cov=fileroot + '/data/proposal.covmat',
#             proposal_scale=1,
#             output_extra_params=['cl_TT']
#         )





#     def __call__(self):
#     	self.cosmo.As = exp(self.cosmo.logA)*1e-10
#         self.cosmo.Yp = self.bbn(**self.cosmo)
#         self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
        
#         self.cmb_result = self.get_cmb(outputs=['cl_TT'],force=True,**self.cosmo)
#         self.cl_TT = {'cl_TT':self.cmb_result['cl_TT'][0:2501]}

#         lnl = self.priors(self)

#         for experiment in self.experiments:
#         	dx = np.dot(experiment.windows, (self.cl_TT)['cl_TT']) - experiment.observation
#         	lnl += -2 * np.dot(dx, np.dot(np.inv(experiment.covariance), dx))

#         return lnl



class gauss_approx(SlikPlugin):
    def __init__(self, *args, **kwargs):
        super(SlikPlugin, self).__init__(*args, **kwargs)
        self.__dict__=self


        try:
            self.observation 
        except:
            self.observation = np.zeros(2501)

        try:
            self.covariance 
        except:
            self.covariance = np.zeros(2501)
	self.inv_cov = np.zeros([2501,2501])
	self.inv_cov[2:2501,2:2501] = inv(self.covariance[2:2501, 2:2501])
        try:
            self.windows
        except:
            self.windows = np.identity(2501)

        try:
            self.cosmo
        except:
            self.cosmo =  get_plugin('models.cosmology')(
                logA = 3.2,
                ns = 0.96,
                ombh2 = 0.0221,
                omch2 = 0.12,
                tau = 0.09,
                theta = 0.0104,
                omnuh2 = 0.000645
                    )
        try:
            self.deltaCl
        except:
            self.deltaCl = np.identity(2501)


        fileroot, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
        self.get_cmb = get_plugin('models.pico')(
            datafile=fileroot + '/data/pico3_tailmonty_v33.dat'
        )

        self.bbn = get_plugin('models.bbn_consistency')()
        self.hubble_theta = get_plugin('models.hubble_theta')() 

    def lnl(self, x):
        self.x_to_cosmo(x)
        self.cosmo['Yp'] = self.bbn(**self.cosmo)
        self.cosmo['H0'] = self.hubble_theta.theta_to_hubble(**self.cosmo)
        self.cosmo.As = np.exp(self.cosmo.logA)*1e-10
        cl_TT = self.get_cmb(outputs=['cl_TT'],force=True,**self.cosmo)['cl_TT'][:2501]
        lnl = self.tau_prior()
        dx = np.dot(self.windows, cl_TT) - self.observation
        lnl += np.dot(dx[2:2501], np.dot(self.inv_cov[2:2501,2:2501], dx[2:2501]))
#	print lnl

        return lnl

    def tau_prior(self):
        dx = self.cosmo.tau - 0.0925
        return np.dot(dx,dx)/0.0125**2

    def x_to_cosmo(self,x):
        self.cosmo['logA'] = x[0]
        self.cosmo['ns'] = x[1]
        self.cosmo['ombh2'] = x[2]
        self.cosmo['omch2'] = x[3]
        self.cosmo['tau'] = x[4]
        self.cosmo['theta'] = x[5]
        self.cosmo['omnuh2'] = x[6]

    def cosmo_to_x(self):
        x = np.zeros(7)
        x[0] = self.cosmo['logA']
        x[1] = self.cosmo['ns'] 
        x[2] = self.cosmo['ombh2'] 
        x[3] = self.cosmo['omch2']
        x[4] = self.cosmo['tau'] 
        x[5] = self.cosmo['theta']
        x[6] = self.cosmo['omnuh2']    
        return x  

    def get_bestfit(self):
        x0 = self.cosmo_to_x()
        minresult = minimize(self.lnl, x0, method = 'Nelder-Mead')
        x_bestfit = minresult.x
	self.x_to_cosmo(x_bestfit)
        self.best_fit = self.get_cmb(outputs=['cl_TT'], force = True, **self.cosmo)['cl_TT'][:2501]

    def sample(self):
        covariance = inv(np.dot(self.deltaCl.T, np.dot(self.inv_cov, self.deltaCl)))
        return self.best_fit + np.dot(self.deltaCl, normal(np.zeros(6), covariance))


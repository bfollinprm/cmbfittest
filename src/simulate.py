import numpy as np
import os.path as osp

from numpy.random import multivariate_normal as normal
from numpy.linalg import inv
from cosmoslik import param_shortcut, get_plugin, SlikDict, SlikPlugin, Slik
from scipy.optimize import minimize



# param = param_shortcut('start','scale')


# class MCMC_chain(SlikPlugin):
#   def __init__(self, experiments):
#         super(SlikPlugin,self).__init__()
#         self.__dict__ = self

#         self.experiments = experiments

#       try: 
#           self.experiments[0]
#       except:
#           self.experiments = [self.experiments]


#         self.cosmo = get_plugin('models.cosmology')(
#             logA = param(3.2),
#             ns = param(0.96),
#             ombh2 = param(0.0221),
#             omch2 = param(0.12),
#             tau = param(0.09,min=0,gaussian_prior=(0.0925,0.015)),
#             theta = param(0.010413),
#             omnuh2 = 0.000645,
#         )
                      


#       fileroot, extension =  osp.split(osp.dirname(osp.abspath(__file__)))


#         self.get_cmb = get_plugin('models.pico')(
#             datafile=fileroot + '/data/pico3_tailmonty_v33.dat'
#         )

#         self.bbn = get_plugin('models.bbn_consistency')()
#         self.hubble_theta = get_plugin('models.hubble_theta')()  
#         self.priors = get_plugin('likelihoods.priors')(self)


#         self.sampler = get_plugin('samplers.metropolis_hastings')(
#             self,
#           num_samples=10000,
#             proposal_cov=fileroot + '/data/proposal.covmat',
#             proposal_scale=1,
#             output_extra_params=['cl_TT']
#         )





#     def __call__(self):
#       self.cosmo.As = exp(self.cosmo.logA)*1e-10
#         self.cosmo.Yp = self.bbn(**self.cosmo)
#         self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
        
#         self.cmb_result = self.get_cmb(outputs=['cl_TT'],force=True,**self.cosmo)
#         self.cl_TT = {'cl_TT':self.cmb_result['cl_TT'][0:self.ellmax + 1]}

#         lnl = self.priors(self)

#         for experiment in self.experiments:
#           dx = np.dot(experiment.windows, (self.cl_TT)['cl_TT']) - experiment.observation
#           lnl += -2 * np.dot(dx, np.dot(np.inv(experiment.covariance), dx))

#         return lnl



class gauss_approx(SlikPlugin):
    def __init__(self, *args, **kwargs):
        super(SlikPlugin, self).__init__(*args, **kwargs)
        self.__dict__=self

        try:
            self.ellmax 
        except:
            self.ellmax = 2500

        try:
            self.observation 
        except:
            self.observation = [np.zeros(self.ellmax + 1)]

        try:
            self.covariance 
        except:
            self.covariance = [np.zeros(self.ellmax + 1)]


        try:
            self.windows
        except:
            self.windows = [np.identity(self.ellmax + 1)]



        try:
            self.deltaCl
        except:
            self.deltaCl = np.identity(self.ellmax + 1)


        self.cosmo =  get_plugin('models.cosmology')(
                logA = 3.2,
                ns = 0.96,
                ombh2 = 0.0221,
                omch2 = 0.12,
                tau = 0.09,
                theta = 0.01046,
                omnuh2 = 0.000645
                )

        self.numexperiments = len(self.covariance)

        self.inv_cov = []
        for i in np.arange(self.numexperiments):
            self.inv_cov.append(inv(self.covariance[i]))

        fileroot, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
        self.get_cmb = get_plugin('models.pico')(
            datafile=fileroot + '/data/pico3_tailmonty_v33.dat'
        )

        self.bbn = get_plugin('models.bbn_consistency')()
        #self.hubble_theta = get_plugin('models.hubble_theta')() 

    def lnl(self, x):
        self.x_to_cosmo(x)
        self.cosmo['Yp'] = self.bbn(**self.cosmo)
        #self.cosmo['H0'] = self.hubble_theta.theta_to_hubble(**self.cosmo)
        self.cosmo.As = np.exp(self.cosmo.logA)*1e-10
        cl_TT = self.get_cmb(outputs=['cl_TT'],force=True,**self.cosmo)['cl_TT'][:self.ellmax + 1]
        lnl = self.tau_prior()
        for i in np.arange(self.numexperiments):
            inv_cov = self.inv_cov[i]
            dx = np.dot(self.windows[i], cl_TT) - self.observation[i]
            if i ==0:
                dx = dx[:-10]
                inv_cov = inv_cov[:-10,:-10]
            temp = np.dot(dx, np.dot(inv_cov, dx))
            # print 'experiment %i gives'%i, temp
            lnl += temp
        # print 'total lnl is', lnl
#   print lnl
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
        self.best_fit = self.get_cmb(outputs=['cl_TT'], force = True, **self.cosmo)['cl_TT'][:self.ellmax + 1]
        #self.best_fit = np.array([self.best_fit for i in np.arange(self.numexperiments)]).flatten()

    def sample(self):
        inv_covariance = np.diag(np.zeros([self.deltaCl.shape[1]]))
        for i in np.arange(self.numexperiments):
            deltaCl = np.dot(self.deltaCl.T, self.windows[i].T)
            inv_covariance += np.dot(deltaCl,np.dot(self.inv_cov[i], deltaCl.T))
        covariance = inv(inv_covariance)
        sample = np.dot(self.deltaCl, normal(np.zeros(6), covariance))
        return sample


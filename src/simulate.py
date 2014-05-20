import numpy as np

from numpy.random import normal

from cosmoslik import param_shortcut, get_plugin, SlikDict, SlikPlugin, Slik



param = param_shortcut('start','scale')


class lcdm_cloud(SlikPlugin):
	def __init__(self, experiments):
        super(SlikPlugin,self).__init__()
        self.__dict__ = self

        self.experiments = experiments

		try: 
			self.experiments[0]
		except:
			self.experiments = [self.experiments]


        self.cosmo = get_plugin('models.cosmology')(
            logA = param(3.2),
            ns = param(0.96),
            ombh2 = param(0.0221),
            omch2 = param(0.12),
            tau = param(0.09,min=0,gaussian_prior=(0.0925,0.015)),
            theta = param(0.010413),
            omnuh2 = 0.000645,
        )
                      


		fileroot, extension =  osp.split(osp.dirname(osp.abspath(__file__)))


        self.get_cmb = get_plugin('models.pico')(
            datafile=fileroot + '/data/pico3_tailmonty_v33.dat'
        )

        self.bbn = get_plugin('models.bbn_consistency')()
        self.hubble_theta = get_plugin('models.hubble_theta')()  
        self.priors = get_plugin('likelihoods.priors')(self)


        self.sampler = get_plugin('samplers.metropolis_hastings')(
            self,
	     	num_samples=10000,
            proposal_cov=fileroot + '/data/proposal.covmat',
            proposal_scale=1,
            output_extra_params=['cl_TT']
        )





    def __call__(self):
    	self.cosmo.As = exp(self.cosmo.logA)*1e-10
        self.cosmo.Yp = self.bbn(**self.cosmo)
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
        
        self.cmb_result = self.get_cmb(outputs=['cl_TT'],force=True,**self.cosmo)
        self.cl_TT = {'cl_TT':self.cmb_result['cl_TT'][0:2501]}

        lnl = self.priors(self)

        for experiment in self.experiments:
        	dx = np.dot(experiment.windows, (self.cl_TT)['cl_TT']) - experiment.observation
        	lnl += -2 * np.dot(dx, np.dot(np.inv(experiment.covariance), dx))

        return lnl




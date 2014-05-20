
class chain(dict):
	def __init__(data, model, *args, **kwargs):
		super(chain,self).__init__(*args, **kwargs)

		try: 
			self.data[0]
		except:
			self.data = [self.data]

		try: 
			self.model[0]
		except:
			self.model = [self.model]

		self.fiducial_cl = (self.model[0]).fiducial_cl



		##Define the cosmology
		modes = []
		for model in self.model:
			modes.append(list(model.dcldphi))
		self.modes = modes

	def lnl(self):
		lnl = 0
		for experiment in self.data:
			theory_bandpowers = dot(experiment.windows, self.cl_TT
			lnl += 2 * dot((experiment.observation - theory_bandpowers).T, dot(inv(experiment.covariance), (experiment.observation - theory_bandpowers)))
		lnl += self.priors()


	def priors(self):

	def cl_TT(self):
		self.cl_TT = self.fiducial_cl + dot(array(self.modes), self.params)
	

	def __call__():	## returns cl_TT
import numpy
import os.path as osp
import compute

class data(dict):
	'''
	class to store an experiment, which is defined as a noise model, window functions, and a set of bandpowers.
	'''

	def __init__(self, experiment, *args, **kwargs):
		super(data, self).__init__(*args, **kwargs)
		self.__dict__ = self
		self.experiment = experiment

	def get_covariance(self):
		root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
		filename = root + '/data/cl_fid.npz'
		try:
			fiducial_model = numpy.load(filename)
			self.fiducial_cl = fiducial_model['cl_TT']
		except:
			raise Exception('Critical Error: No fiducial model found, or fiducial model in incorrect format.  Need file %s to contain dictionary with keys cl_TT and params'%filename)

		self.covariance = compute.calc_covariance(self)


	def get_window_functions(self):
		self.windows = compute.load_window_functions(self)


	def get_bandpowers(self):
		try:
			self.observation
		except:

			self.observation = compute.get_bandpowers(self)



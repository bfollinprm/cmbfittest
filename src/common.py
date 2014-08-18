
#import simulate





# def calc_covariance(data):


# 	root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
# 	fileroot = root + '/data/experiments/%s/'%data.experiment

# 	if data.experiment == 'spt':
# 		filename = fileroot + 'covariance.dat'
# 		covariance = np.loadtxt(filename)
		

# 	if data.experiment == 'wmap':
# 		fidCl = data.fiducial_cl[:data.ellmax + 1]
# 		ondiag = np.loadtxt(fileroot + '/data/highl/wmap_likelihood_inputs_tt.p4v6.wmap9.kq85.cinv_v3.dat')
# 		offdiag = np.loadtxt(fileroot + '/data/highl/wmap_likelihood_inputs_tt_offdiag.p4v4.wmap9.kq85.cinv_v3.dat')


# 		ell = ondiag[:,0].astype('int')
# 		bandpower = ondiag[:,1]
# 		Nl = ondiag[:,2]
# 		Fsky = ondiag[:,3]
# 		l1 = offdiag[:,0].reshape([1199, 599]).astype('int')
# 		l2 = offdiag[:,1].reshape([1199, 599]).astype('int')
# 		Epsilon_temp = offdiag[:,2].reshape([1199, 599])
# 		R_eff_temp = offdiag[:,3].reshape([1199, 599])

# 		#R_eff = np.diag(np.zeros(data.ellmax + 1))
# 		R_eff = np.diag(np.zeros(1201))
# 		R_eff[l1, l2] = R_eff_temp
# 		R_eff[l2, l1] = R_eff_temp


# 		Epsilon = np.diag(np.zeros(1201))
# 		#np.diag(np.zeros(data.ellmax + 1))
# 		Epsilon[l1, l2] = Epsilon_temp
# 		Epsilon[l2, l1] = Epsilon_temp
# 		tttt = np.zeros(1201)
# 		#np.zeros(data.ellmax + 1) + 1
	
# 		tttt[ell] = 2 * (fidCl[ell] + Nl)**2/(2*ell + 1)/Fsky**2


# 		tttt = np.outer(np.sqrt(tttt), np.sqrt(tttt))

# 		inv_covariance = np.diag(1.0/np.diag(tttt)) - R_eff/(tttt) + Epsilon/tttt**2
# 		covariance = np.diag(np.zeros(data.ellmax + 1)+1.0e30)
# 		covariance[2:1201, 2:1201] = inv(inv_covariance[2:,2:])


	# if data.experiment == 'planck_like':
	# 	ell = np.arange(data.ellmax + 1)
	# 	Fsky = 0.3
	# 	Nl  = 0.00015048 * np.exp(ell**2 * 7.69112e-7 ) * (ell**2) / 2.0 / np.pi
	# 	covariance = 2 * (data.fiducial_cl[:data.ellmax + 1] + Nl)**2/(2*ell + 1)/Fsky
	# 	covariance = np.diag(covariance)
	# return covariance



# def load_window_functions(data):

# 	root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
# 	fileroot = root + '/data/experiments/%s/'%data.experiment

# 	if data.experiment == 'spt':

# 		windows = [np.loadtxt(fileroot + 'windows/window_%i'%i) for i in range(1,48)]
# 		ell = np.array(windows)[0,:,0].astype('int')
# 		weights = np.array(windows)[:,:,1]
# 		windows = np.zeros([47,data.ellmax + 1])
# 		ell = [i for i in ell if i < data.ellmax + 1]
# 		windows[:, ell] = weights[:,np.array(ell)-50]
# 		return windows

# 	if data.experiment == 'wmap':

# 		windows = np.identity(data.ellmax + 1)[:1200,:]
# 		return windows

# 	if data.experiment == 'planck_like':

# 		windows = np.identity(data.ellmax + 1)
# 		return windows

# 	return windows


# def get_bandpowers(data):
# 	root, extension =  osp.split(osp.dirname(osp.abspath(__file__)))
# 	fileroot = root + '/data/experiments/%s/'%data.experiment

# 	if data.experiment == 'spt':

# 		bandpowers = np.loadtxt(fileroot + 'bandpowers.dat')

# 	if data.experiment == 'wmap':

# 		ondiag = fileroot + 'wmap_likelihood_inputs_tt.p4v6.wmap9.kq85.cinv_v3.dat'
# 		ell = ondiag[:,0]
# 		bandpowers_tmp = ondiag[:,1]
# 		bandpowers = np.zeros(data.ellmax + 1)
# 		bandpowers[ell] = bandpowers_tmp


# 	if data.experiment == 'planck_like':

# 		bandpowers = data.fiducial_cl[:data.ellmax + 1] + normal(np.zeros(data.covariance.shape[0]), data.covariance)

# 	return bandpowers


# def get_modes(test):

# 	if test.type in ['LDA', 'PCA', 'S_NL', 'S_N']:
# 		test_matrix = np.dot(inv(test.noise_cov[2:,2:] + test.kappa * test.SM_cov[2:,2:]), test.supermodel_cov[2:,2:])
# 		w, modes = eig(test_matrix)
	
# 	if test.type in ['fisher', 'FN', 'FNL']:
# 		test_matrix = np.dot(inv(test.noise_cov[2:,2:]) + test.kappa + test.SM_cov[2:,2:])
# 		w, modes = eig(test_matrix)

# 	padded_modes = np.zeros([test.ellmax + 1, test.ellmax - 1])
# 	for i in np.arange((np.real(modes).shape)[0]):
# 		padded_modes[2:,i] = np.real(modes[:,i])

# 	##Make sure they're sorted
# 	#ii = np.argsort(w)
# 	#w = w[ii]
# 	#w = w[::-1]
# 	#padded_modes = modes[:,ii]
# 	#padded_modes = modes[:,::-1]
# 	return np.real(w), np.real(padded_modes) 

	




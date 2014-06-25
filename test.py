import cmbfittest as fit
from pylab import *
#import numba


from scipy.ndimage.filters import gaussian_filter1d as gaussian_filter
from numpy.random import multivariate_normal as normal
from numpy.random import normal as N
from numpy import dot, sqrt, zeros, array, arange, save, mean, std
from numpy.linalg import inv, eigh
from scipy.linalg import sqrtm


##Instance model structures for linearized LCDM and primordial perturbations
print 'instancing models'
primordial = fit.model.model('primordial')
lcdm = fit.model.model('lcdm')


##Instance experimental noise model
experiment = fit.data.data('planck_like')
experiment.get_covariance()

##Calculate or load linear responses to model parameters
primordial.get_dcldphi()
lcdm.get_dcldphi()


###We don't care about ell > 2500, truncate for speed-up
primordial.dcldphi = primordial.dcldphi[:2501,:]
lcdm.dcldphi = lcdm.dcldphi[:2501,:]


###Numerics mean the responses in Cl space to perturbations are rocky; smooth them out
lcdm.dcldphi = gaussian_filter(lcdm.dcldphi, sigma = 20, axis = 0)
primordial.dcldphi= gaussian_filter(primordial.dcldphi, sigma = 20, axis = 0)


###One choice of prior for LCDM modes--everything has the same weight
print '\t ...setting priors'
unit_lcdm = fit.model.model('lcdm')
row_sums = sqrt((lcdm.dcldphi**2).sum(axis = 0))
unit_lcdm.dcldphi = lcdm.dcldphi/row_sums


###Another choice of prior--LCDM variation looks as much like noise as possible (Fisher prior)
fisher_lcdm = fit.model.model('lcdm')
parameter_prior_cov = inv(dot(lcdm.dcldphi[2:,:].T, dot(inv(experiment.covariance[2:2501,2:2501]), lcdm.dcldphi[2:,:])))
fisher_lcdm.dcldphi = dot(lcdm.dcldphi, parameter_prior_cov)




###Get the covariance
print '\t ...calculating covariances'
primordial.get_covariance()
lcdm.get_covariance()
unit_lcdm.get_covariance()
fisher_lcdm.get_covariance()

##Instance a compression scheme
print 'instancing compression scheme'
S_N = fit.scheme.compression_scheme(
                                    'fisher', 
                                    dcldphi = primordial.dcldphi, 
                                    noise_cov = experiment.covariance, 
                                    #supermodel_cov = primordial.covariance,
                                    #kappa = 1000, 
                                    #SM_cov = fisher_lcdm.covariance,
                                    nummodes = 20
                                    )
print '\t ...calculating modes'
S_N.get_compression_modes()


plot(S_N.modes[:,:5])
show()


n = 30
chisqlist = zeros(n)
maxchisqlist = zeros(n)
globalz = zeros(n)
maxz = zeros(n)
for ii in arange(n):

    observation = experiment.fiducial_cl[:2501].copy()
    observation += normal(zeros(2501), experiment.covariance[:2501,:2501])

    print 'simulating spectra %i'%ii
    cloud = fit.sim.gauss_approx(observation = observation, 
                                  covariance = experiment.covariance, 
                                  deltaCl = lcdm.dcldphi)
    print '\t ...finding best fit'
    cloud.get_bestfit()


    ##And sample the theory modes
    print '\t ...sampling'
    draws = array([cloud.sample() for i in arange(10000)])


    ##And add noise (diagonalize noise model first and this goes MUCH faster).

    variance, orthogonal_modes = eigh(experiment.covariance[:2501,:2501])
    sigma = sqrt(variance)

    print '\t ...adding noise'


    #@jit
    def randN(scale = 0):
        if scale == 0:
          return 0
        else:
          return N(scale = scale)

    def add_noise(spectrum, sigma, orthogonal_modes):
      return spectrum + dot(orthogonal_modes, array([randN(scale = sigma[j]) for j in arange(2501)]))

    noisydraws = array(map(lambda x:add_noise(x, sigma, orthogonal_modes), draws.tolist()))

    ###Turn the draws and observation into test amplitudes.

    print 'calculating test statistics'

    cloud_amps, observed_amps = fit.test.get_amplitudes(noisydraws, cloud.observation, S_N.modes)


    chisq, zscore = fit.test.global_statistic(cloud_amps, observed_amps)
    chisqlist[ii] = chisq
    globalz[ii] = zscore
    print '\t ...reduced chi2 is %g'%chisq
    print '\t \t for a z-score of %g'%zscore

    maxchisq, zscore = fit.test.max_statistic(cloud_amps, observed_amps)
    maxchisqlist[ii] = maxchisq
    maxz[ii] = zscore
    print '\t ...max statistic is %g'%maxchisq
    print '\t \t for a z-score of %g'%zscore


expectation = mean(globalz)
error = std(globalz)
print 'Global zscore:  %g'%expectation, ' \pm %g'%error

expectation = mean(maxz)
error = std(maxz)
print 'Max zscore:  %g'%expectation,  ' \pm %g'%error


save('maxchisq', maxchisqlist)
save('chisq', chisqlist)
save('globalz', globalz)
save('maxz', maxz)



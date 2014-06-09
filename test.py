# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import model

# <codecell>

primordial = model.model('primordial')
lcdm = model.model('lcdm')

# <codecell>

primordial.get_dcldphi()
lcdm.get_dcldphi()

# <codecell>

primordial.dcldphi = primordial.dcldphi[:2501,:]
lcdm.dcldphi = lcdm.dcldphi[:2501,:]

# <codecell>

from scipy.ndimage.filters import gaussian_filter1d as gaussian_filter

# <codecell>

lcdm.dcldphi_smooth = gaussian_filter(lcdm.dcldphi, sigma = 20, axis = 0)
primordial.dcldphi_smooth = gaussian_filter(primordial.dcldphi, sigma = 20, axis = 0)

# <codecell>

plot(primordial.dcldphi_smooth);
primordial.save_dcldphi()
lcdm.save_dcldphi() 
#.dcldphi / dot(lcdm.dcldphi, lcdm.dcldphi).reshape(-1,1)

# <codecell>

row_sums = sqrt((lcdm.dcldphi**2).sum(axis = 0))
lcdm.dcldphi = lcdm.dcldphi /row_sums
plot(lcdm.dcldphi);

# <codecell>

primordial.get_covariance()
lcdm.get_covariance()

# <codecell>

imshow(lcdm.covariance)
colorbar()
import data

# <codecell>

experiment = data.data('planck_like')

# <codecell>

experiment.get_covariance()

# <codecell>

import compression_scheme

# <codecell>

lda = compression_scheme.compression_scheme(
                                            'LDA', 
                                            dcldphi = primordial.dcldphi_smooth, 
                                            noise_cov = experiment.covariance, 
                                            supermodel_cov = primordial.covariance 
                                            #kappa = 1.0e10, 
                                            #SM_cov = lcdm.covariance
                                            )

# <codecell>

lda.get_compression_modes()

# <codecell>

plot(lda.modes[5:,19]);

# <codecell>

import simulate
from numpy.random import multivariate_normal as normal
#observation = experiment.fiducial_cl[:2501] + normal(np.zeros(2501), experiment.covariance[:2501,:2501])

# <codecell>

clobs = loadtxt('../clobs.csv', delimiter = ',')
#plot(clobs[:,2:])
observation = experiment.fiducial_cl[:2501].copy()
observation += normal(np.zeros(2501), experiment.covariance[:2501,:2501])
observation2 = experiment.fiducial_cl[:2501].copy()
observation2[clobs[:,0].astype('int')] = clobs[:,4]
plot(observation)

# <codecell>

cloud = simulate.gauss_approx(observation = observation, 
                              covariance = experiment.covariance, 
                              deltaCl = lcdm.dcldphi_smooth)

# <codecell>

cloud.get_bestfit()
plot(cloud.best_fit - experiment.fiducial_cl[:2501])
ylim(-100,100)

# <codecell>

plot(cloud.best_fit)
plot(experiment.fiducial_cl[:2501], alpha = 1)
plot(cloud.observation, alpha = 0.3, color = 'k')


# <codecell>

x = array([cloud.sample() for i in arange(10000)])
x.shape

# <codecell>

#plot(x.T, color = 'b', alpha = .3);
#plot(cloud.observation, color = 'k', alpha = 0.3);
#ylim(1,1e4);

# <codecell>

#plot(array([x[i,:]/cloud.best_fit for i in arange(10000)]).T, color = 'b', alpha = 0.3);
#plot(cloud.observation/cloud.best_fit, color = 'k', alpha = 0.3);
#ylim(0,2);

# <codecell>

#lda.modes.shape
#x.shape
sigma2 = diag(cloud.covariance)
sigma2[0] = 1
xwithnoise = array([x[j,:] + array([random.normal(scale = sqrt(sigma2[i])) for i in arange(0,2501)]) for j in arange(10000)])

# <codecell>

cloud.amps = dot(lda.modes.T, xwithnoise.T)
#cloud.amps.shape

# <codecell>

observed_amp = dot(lda.modes.T, cloud.observation)
#observed_amp.shape

# <codecell>

from scipy.linalg import sqrtm
means = mean(cloud.amps, axis = 1)
covariance = cov(cloud.amps)
covariance.shape

# <codecell>

z = dot(sqrtm(covariance), observed_amp)
zcloud = dot(sqrtm((covariance)), cloud.amps)
zcloud.shape

# <codecell>


# <codecell>

klist = [randint(0,10000) for i in arange(1000)]
distance = [dot((observed_amp[:i] - means[:i]).T, dot(inv(covariance[:i,:i]), (observed_amp[:i] - means[:i]).T)) for i in arange(1,observed_amp.size)]
test_distances = [[dot((cloud.amps[:i,k] - means[:i]).T, dot(inv(covariance[:i,:i]), (cloud.amps[:i,k] - means[:i]).T)) for i in arange(1,observed_amp.size)] for k in klist]

# <codecell>

dis = [(observed_amp[i] - means[i])**2/diag(covariance)[i] for i in arange(1, observed_amp.size)]
test_dis = [[(cloud.amps[i,k] - means[i])**2/diag(covariance)[i] for i in arange(1, observed_amp.size)] for k in arange(1000)]
average = mean(array(test_dis), axis = 0)
sigma = std(array(test_dis),axis = 0)

#errorbar(arange(29), average, yerr= sigma)
plot(array([zip(arange(29), array(test_dis)[k,:]) for k in arange(1000)])[:,:,0].T, array([zip(arange(29), array(test_dis)[k,:]) for k in arange(1000)])[:,:,1].T, 
     color = 'k', alpha = 0.01);
scatter(arange(29), dis, color = 'red')
ylim(0,10)



# <codecell>


# <codecell>

plot(array([dot(lda.modes, cloud.amps[:,i]) - dot(lda.modes, observed_amp) for i in arange(10000)]).T, color = 'k', alpha = .01);
#plot(dot(lda.modes, observed_amp))
ylabel(r'$C_l^{cloud}  - C_l^{observed}$--with noise', fontsize = 20);
title('S_N projected power; no exotic perturbations');

# <codecell>

semilogy(lda.eigspectrum[:30]/max(lda.eigspectrum))

# <codecell>



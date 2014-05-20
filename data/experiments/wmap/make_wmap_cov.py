from pylab import *


ondiag = genfromtxt('data/highl/wmap_likelihood_inputs_tt.p4v6.wmap9.kq85.cinv_v3.dat')
offdiag = genfromtxt('data/highl/wmap_likelihood_inputs_tt_offdiag.p4v4.wmap9.kq85.cinv_v3.dat')


ell = ondiag[:,0]
bandpower = ondiag[:,1]
Nl = ondiag[:,2]
Fsky = ondiag[:,3]

l1 = offdiag[:,0].reshape([1199, 599]).astype('int')
l2 = offdiag[:,1].reshape([1199, 599]).astype('int')
Epsilon_temp = offdiag[:,2].reshape([1199, 599])
R_eff_temp = offdiag[:,3].reshape([1199, 599])

R_eff = diag(zeros(1200))
R_eff[l1-2, l2-2] = R_eff_temp
R_eff[l2-2, l1-2] = R_eff_temp


Epsilon = diag(zeros(1200))
Epsilon[l1-2, l2-2] = Epsilon_temp
Epsilon[l2-2, l1-2] = Epsilon_temp

tttt = 2 * (cl_fid[ell] + Nl)**2/(2*ell + 1)/Fsky**2

tttt = numpy.outer(sqrt(tttt), sqrt(tttt))


Fisher = diag(1/diag(tttt)) - R_eff/sqrt(tttt) + Epsilon/tttt
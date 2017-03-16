import numpy as np

import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla

N = 5812
Re = 80

savedmatsstr = '../data/cylinderwake__mats_NV{1}_Re{0}.mat'.\
    format(1, N)
print 'Re = ', Re
print 'NV = ', N

# get the coefficients
mats = scipy.io.loadmat(savedmatsstr)
A = 1./Re*mats['A']
M = mats['M']
J = mats['J']
fv = mats['fv'] + 1./Re*mats['fv_diff']
fp = mats['fp'] + mats['fp_div']
NV, NP = fv.shape[0], fp.shape[0]

stksrhs = np.vstack([fv, fp])
stksmat = sps.vstack([sps.hstack([A, -J.T]),
                      sps.hstack([J, sps.csc_matrix((NP, NP))])])
stksvp = spsla.spsolve(stksmat, stksrhs).reshape((NV+NP, 1))
stksv = stksvp[:NV].reshape((NV, 1))
stksp = stksvp[NV:].reshape((NP, 1))

resstksmom = A*stksv - J.T*stksp - fv
resconti = J*stksv - fp
print 'The Stokes momentum eq residual: ', np.linalg.norm(resstksmom)
print 'The conti residual: ', np.linalg.norm(resconti)

try:
    savedmatsstr = '../data/cylinderwake__mats_NV{1}_Re{0}.mat'.\
        format(1, N)
    mats = scipy.io.loadmat(savedmatsstr)
    vssstks = mats['v_ss_stokes']
    pssstks = mats['p_ss_stokes']
    print 'Difference of v solutions: ', np.linalg.norm(stksv - vssstks)
    print 'Difference of p solutions: ', np.linalg.norm(stksp - pssstks)
except IOError:
    print 'No data for Re={0}, N={1}'.format(Re, N)

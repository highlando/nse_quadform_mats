import numpy as np

import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla

import conv_tensor_utils as ctu

N = 2
Re = 40
palpha = 1e-3  # penalization parameter for the robin relaxation
uvec = 1*np.array([[1], [-1]])  # the steady input

NVdict = {1: 5812, 2: 9384}

NV = NVdict[N]
savedmatsstr = '../data/cylinderwake__mats_' +\
    'NV{1}_Re{0}_bccontrol_palpha{2}.mat'.format(1, NV, 1)
visujsonstr = '../data/visualization_cylinderwake_N{0}.jsn'.format(N)

print 'Re = ', Re
print 'NV = ', NV

npicardstps = 5

# get the coefficients
mats = scipy.io.loadmat(savedmatsstr)
A = 1./Re*mats['A'] + mats['L1'] + mats['L2'] + 1./palpha*mats['Arob']
M = mats['M']
J = mats['J']
hmat = mats['H']
Brob = mats['Brob']

fv = mats['fv'] + 1./Re*mats['fv_diff'] + mats['fv_conv'] \
    + 1./palpha*np.dot(Brob, uvec)
fp = mats['fp'] + mats['fp_div']
NV, NP = fv.shape[0], fp.shape[0]

updnorm = 1
curv = np.zeros((NV, 1))
curvp = np.zeros((NV+NP, 1))

stpcount = 0
while updnorm > 1e-10:
    picard = True if stpcount < npicardstps else False

    def linnsevp(nxvp):
        nxv = nxvp[:NV].reshape((NV, 1))
        nxp = nxvp[NV:].reshape((NP, 1))
        if picard:
            momeq = A*nxv + hmat*np.kron(nxv, curv) - J.T*nxp
        else:
            momeq = A*nxv + hmat*np.kron(nxv, curv) + hmat*np.kron(curv, nxv) \
                - J.T*nxp
        contieq = J*nxv
        return np.vstack([momeq.reshape((NV, 1)),
                          contieq.reshape((NP, 1))])

    aph = spsla.LinearOperator((NV+NP, NV+NP), matvec=linnsevp)

    if picard:
        currhs = np.vstack([fv, fp])
    else:
        currhs = np.vstack([fv+hmat*np.kron(curv, curv), fp])

    H1k, H2k = ctu.linearzd_quadterm(hmat, curv, retparts=True)
    HL = H1k if picard else H1k+H2k
    cursysmat = sps.vstack([sps.hstack([A+HL, -J.T]),
                            sps.hstack([J, sps.csc_matrix((NP, NP))])])
    nextvp = spsla.spsolve(cursysmat, currhs).reshape((NV+NP, 1))

    ittype = 'Picard' if picard else 'Newton'
    print 'Iteration step {0} ({1})'.format(stpcount, ittype)
    nextv = nextvp[:NV].reshape((NV, 1))
    nextp = nextvp[NV:].reshape((NP, 1))
    curnseres = A*nextv + hmat*np.kron(nextv, nextv) - J.T*nextp - fv
    print 'norm of nse residual', np.linalg.norm(curnseres)
    updnorm = np.linalg.norm(nextv - curv) / np.linalg.norm(nextv)
    print 'norm of current update :', updnorm
    curv = nextv
    stpcount += 1

print '\n *** Done *** \n'
print 'NSE momentum eq residual: ', np.linalg.norm(curnseres)
resconti = J*nextv - fp
print 'The conti residual: ', np.linalg.norm(resconti)

import visualization_utils as vu
vu.cylwake_paraview(pvec=nextp, velvec=nextv, strtojson=visujsonstr,
                    pfile='p__', vfile='vel__')

import numpy as np
import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import conv_tensor_utils as ctu
import visualization_utils as vu

# hard coded paths and dictory for data
NVdict          = {1: 5824, 2: 9384,  3: 19512}
savedmatsstr    = lambda NV: '../data/cylinderwake__mats_NV{1}_Re{0}_bccontrol_palpha{2}.mat'.format(1,NV,1)
visujsonstr     = lambda N : '../data/visualization_cylinderwake_N{0}.jsn'.format(N)


# setup parameters
Re          = 40
palpha      = 1e-3                      # penalty for robin boundary control
uvec        = 1*np.array([[1], [-1]])   # the steady input
npicardstps = 5
N           = 1
NV          = NVdict[N]


# visualisation files
pfile = 'p__cylinderwake_stst_bccontrol_Re{0}_NV{1}_palpha{2:e}.vtu'.format(Re, NV, palpha)
vfile = 'v__cylinderwake_stst_bccontrol_Re{0}_NV{1}_palpha{2:e}.vtu'.format(Re, NV, palpha)


# print reynolds number and discretization lvl
print 'Re = {0}'.format(Re)
print 'NV = {0}'.format(NV)
print '\n'

# load the coefficients matrices
mats    = scipy.io.loadmat(savedmatsstr(NV))
M       = mats['M']
A       = 1./Re*mats['A'] + mats['L1'] + mats['L2'] + 1./palpha*mats['Arob']
J       = mats['J']
hmat    = mats['H']
Brob    = mats['Brob']
fv      = mats['fv'] + 1./Re*mats['fv_diff'] + mats['fv_conv'] + 1./palpha*np.dot(Brob,uvec)
fp      = mats['fp'] + mats['fp_div']
NV, NP  = fv.shape[0], fp.shape[0]


# solve steady state equation
updnorm     = 1
curv        = np.zeros((NV, 1))
curvp       = np.zeros((NV+NP, 1))
stpcount    = 0

while updnorm > 1e-10:
    picard = True if stpcount < npicardstps else False

    if picard:
        currhs = np.vstack([fv, fp])
    else:
        currhs = np.vstack([fv+hmat*np.kron(curv, curv), fp])

    H1k, H2k    = ctu.linearzd_quadterm(hmat, curv, retparts=True)
    HL          = H1k if picard else H1k+H2k
    cursysmat   = sps.vstack([sps.hstack([A+HL, -J.T]),sps.hstack([J, sps.csc_matrix((NP, NP))])]).tocsc()
    nextvp      = spsla.spsolve(cursysmat, currhs).reshape((NV+NP, 1))

    print 'Iteration step {0} ({1})'.format(stpcount, 'Picard' if picard else 'Newton')
    nextv       = nextvp[:NV].reshape((NV, 1))
    nextp       = nextvp[NV:].reshape((NP, 1))
    curnseres   = A*nextv + hmat*np.kron(nextv, nextv) - J.T*nextp - fv
    print 'Norm of nse residual:   {0:e}'.format(np.linalg.norm(curnseres))
    updnorm     = np.linalg.norm(nextv - curv) / np.linalg.norm(nextv)
    print 'Norm of current update: {0:e}'.format(updnorm)
    print '\n'
    curv = nextv
    stpcount += 1


# print results
print '*** Done ***'
print 'NSE momentum eq residual: {0:e}'.format(np.linalg.norm(curnseres))
resconti = J*nextv - fp
print 'The conti residual:       {0:e}'.format(np.linalg.norm(resconti))
print '\n'

# write paraview
vu.writevp_paraview(pvec=nextp, velvec=nextv, strtojson=visujsonstr(N),pfile=pfile, vfile=vfile)
print '*** for visualization try ***'
print 'paraview {0}'.format(vfile)
print 'paraview {0}'.format(pfile)




import numpy as np
import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import conv_tensor_utils as ctu
import visualization_utils  as vu
import sys, getopt

# hard coded paths and dictionary for data
NVdict          = {10: 722, 20: 3042,  30: 6962}
savedmatsstr    = lambda NV: '../data/drivencavity__mats_NV{1}_Re{0}.mat'.format(1,NV)
visujsonstr     = lambda NV : '../data/visualization_drivencavity_NV{0}.jsn'.format(NV)


# setup parameters
N           = 30
Re          = 40
npicardstps = 5


# get command line input and overwrite standard paramters if necessary
options, rest = getopt.getopt(sys.argv[1:], '',['N=', 'Re=', 'Picardsteps='])
for opt, arg in options: 
    if opt == '--N':
        N = int(arg)
    elif opt == '--Re':
        Re = int(arg)
    elif opt == '--Picardsteps':
        npicardstps = int(arg)


# visualisation files
NV    = NVdict[N]
pfile = 'results/p__drivencavity_stst_Re{0}_NV{1}.vtu'.format(Re, NV)
vfile = 'results/v__drivencavity_stst_Re{0}_NV{1}.vtu'.format(Re, NV)


# print reynolds number and discretization lvl
print('Re           = {0}'.format(Re))
print('NV           = {0}'.format(NV))
print('Picardsteps  = {0}'.format(npicardstps))
print('pfile        = {0}'.format(pfile))
print('vfile        = {0}'.format(vfile))
print('\n')


# load the coefficients matrices
mats    = scipy.io.loadmat(savedmatsstr(NV))
M       = mats['M']
A       = 1./Re*mats['A'] + mats['L1'] + mats['L2']
J       = mats['J']
hmat    = mats['H']
fv      = mats['fv'] + 1./Re*mats['fv_diff'] + mats['fv_conv']
fp      = mats['fp'] + mats['fp_div']
NV, NP  = fv.shape[0], fp.shape[0]


# solve steady state equation
updnorm     = 1
curv        = np.zeros((NV, 1))
curvp       = np.zeros((NV+NP, 1))
stpcount    = 0

while updnorm > 1e-10:

    H1k, H2k    = ctu.linearzd_quadterm(hmat, curv, retparts=True)
    picard      = stpcount < npicardstps
    if picard:
        currhs = np.vstack([fv, fp])
        HL     = H1k
    else:
        currhs = np.vstack([fv+ctu.eva_quadterm(hmat,curv), fp])
        HL     = H1k + H2k

    cursysmat   = sps.vstack([sps.hstack([A+HL, -J.T]),sps.hstack([J, sps.csc_matrix((NP, NP))])]).tocsc()
    nextvp      = spsla.spsolve(cursysmat, currhs).reshape((NV+NP, 1))

    print('Iteration step {0} ({1})'.format(stpcount, 'Picard' if picard else 'Newton'))
    nextv       = nextvp[:NV].reshape((NV, 1))
    nextp       = nextvp[NV:].reshape((NP, 1))
    curnseres   = A*nextv + ctu.eva_quadterm(hmat,nextv) - J.T*nextp - fv
    print('Norm of nse residual:   {0:e}'.format(np.linalg.norm(curnseres)))
    updnorm     = np.linalg.norm(nextv - curv) / np.linalg.norm(nextv)
    print('Norm of current update: {0:e}'.format(updnorm))
    print('\n')
    curv = nextv
    stpcount += 1


# print results
print('*** Done ***')
print('NSE momentum eq residual: {0:e}'.format(np.linalg.norm(curnseres)))
resconti = J*nextv - fp
print('The conti residual:       {0:e}'.format(np.linalg.norm(resconti)))
print('\n')


# write paraview
vu.writevp_paraview(pvec=nextp, velvec=nextv, strtojson=visujsonstr(NV),pfile=pfile, vfile=vfile)
print('*** for visualization try ***')
print('paraview {0}'.format(vfile))
print('paraview {0}'.format(pfile))





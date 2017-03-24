import numpy as np
import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import sys, getopt


# hard coded paths and dictionary for data
NVdict          = {1: 5812, 2: 9356,  3: 19468}
savedmatsstr    = lambda NV: '../data/cylinderwake__mats_NV{1}_Re{0}.mat'.format(1,NV)
visujsonstr     = lambda NV : '../data/visualization_cylinderwake_NV{0}.jsn'.format(NV)


# setup parameters
N           = 1
Re          = 80
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


# futher parameters
NV          = NVdict[N]

# visualisation files
pfile = 'results/p__cylinderwake_stst_Re{0}_NV{1}.vtu'.format(Re, NV)
vfile = 'results/v__cylinderwake_stst_Re{0}_NV{1}.vtu'.format(Re, NV)


#create dir if not exists
if not os.path.exists('results'):
    os.makedirs('results')


# print reynolds number and discretization lvl
print('Re           = {0}'.format(Re))
print('NV           = {0}'.format(NV))
print('Picardsteps  = {0}'.format(npicardstps))
print('\n')


# load the coefficients matrices
mats    = scipy.io.loadmat(savedmatsstr(NV))
M       = mats['M']
A       = 1./Re*mats['A']
J       = mats['J']
hmat    = mats['H']
fv      = mats['fv'] + 1./Re*mats['fv_diff']
fp      = mats['fp'] + mats['fp_div']
NV, NP  = fv.shape[0], fp.shape[0]


# solve stokes equation
stksrhs = np.vstack([fv, fp])
stksmat = sps.vstack([sps.hstack([A, -J.T]),sps.hstack([J, sps.csc_matrix((NP, NP))])]).tocsc()
stksvp  = spsla.spsolve(stksmat, stksrhs).reshape((NV+NP, 1))
stksv   = stksvp[:NV].reshape((NV, 1))
stksp   = stksvp[NV:].reshape((NP, 1))


# print results
resstksmom  = A*stksv - J.T*stksp - fv
resconti    = J*stksv - fp
print('*** Done ***')
print('The Stokes momentum eq residual: {0:e}'.format(np.linalg.norm(resstksmom)))
print('The conti residual:              {0:e}'.format(np.linalg.norm(resconti)))

vssstks = mats['v_ss_stokes']
pssstks = mats['p_ss_stokes']
 
print('Difference of v solutions:       {0:e}'.format(np.linalg.norm(stksv - vssstks)))
print('Difference of p solutions:       {0:e}'.format(np.linalg.norm(stksp - pssstks)))


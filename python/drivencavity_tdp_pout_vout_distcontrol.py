import numpy as np

import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla

import conv_tensor_utils as ctu
import visualization_utils as vu

N = 20
Re = 800

NVdict = {10: 722, 20: 3042, 30: 6962}

NV = NVdict[N]
savedmatsstr = '../data/drivencavity__mats_NV{1}_Re{0}.mat'.\
    format(1, NV)
visujsonstr = '../data/visualization_drivencavity_NV{0}.jsn'.format(NV)
print 'Re = ', Re
print 'NV = ', NV

# get the coefficients
mats = scipy.io.loadmat(savedmatsstr)
A = 1./Re*mats['A'] + mats['L1'] + mats['L2']
B = mats['B']
M = mats['M']
J = mats['J']
hmat = mats['H']
fv = mats['fv'] + 1./Re*mats['fv_diff'] + mats['fv_conv']
fp = mats['fp'] + mats['fp_div']
pcmat = mats['Cp']
vcmat = mats['Cv']
NV, NP = fv.shape[0], fp.shape[0]

# Fix the p
J = J[:-1, :]
fp = fp[:-1, :]

t0, tE, Nts = 0., 20., 2**12
trange = np.linspace(t0, tE, Nts+1)
DT = (tE-t0)/Nts
omeg = 4.  # parameter for the frequency of the input signal

rdir = 'results/'
vfileprfx = 'v_tdpdrivcav_NV{0}_Re{1}_tE{2}_Nts{3}_bccomg{4}'.\
    format(NV, Re, tE, Nts, omeg)
pfileprfx = 'p_tdpdrivcav_NV{0}_Re{1}_tE{2}_Nts{3}_bccomg{4}'.\
    format(NV, Re, tE, Nts, omeg)

# restrict it to two dofs in the input
NU = B.shape[1]
B = B[:, [0, NU/2]]


# define an input function
def bbcu(t):
    uvec = np.array([[np.sin(omeg*np.pi*t/(tE-t0))],
                    [np.cos(omeg*np.pi*t/(tE-t0))]])
    return B*uvec

sysmat = sps.vstack([sps.hstack([M+DT*A, -J.T]),
                     sps.hstack([J, sps.csc_matrix((NP-1, NP-1))])])

print 'computing Stokes solution to be used as initial value...'

fvstks = mats['fv'] + 1./Re*mats['fv_diff'] + bbcu(t0)
Astks = 1./Re*mats['A']

stksrhs = np.vstack([fvstks, fp])
stksmat = sps.vstack([sps.hstack([Astks, -J.T]),
                      sps.hstack([J, sps.csc_matrix((NP-1, NP-1))])])
stksvp = spsla.spsolve(stksmat, stksrhs).reshape((NV+NP-1, 1))
stksv = stksvp[:NV].reshape((NV, 1))
stksp = (np.r_[stksvp[NV:].flatten(), 0]).reshape((NP, 1))

print 'computing LU once...'
sysmati = spsla.factorized(sysmat)

print 'doing the time loop...'
old_v = stksv

# Preparing for the output
poutlist = []
voutlist = []
vfile = rdir + vfileprfx + '__t{0}.vtu'.format(trange[0])
pfile = rdir + pfileprfx + '__t{0}.vtu'.format(trange[0])
vu.writevp_paraview(velvec=stksv, pvec=stksp, vfile=vfile, pfile=pfile,
                    strtojson=visujsonstr)
vfilelist, pfilelist = [vfile], [pfile]

for k, t in enumerate(trange):
    crhsv = M*old_v + DT*(fv - ctu.eva_quadterm(hmat, old_v) + bbcu(t))
    crhs = np.vstack([crhsv, fp])
    vp_new = np.atleast_2d(sysmati(crhs.flatten())).T
    old_v = vp_new[:NV]
    # p = vp_new[NV:]
    p = (np.r_[vp_new[NV:].flatten(), 0]).reshape((NP, 1))

    poutlist.append((pcmat*p)[0][0])
    voutlist.append((vcmat*old_v).flatten())
    if np.mod(k, Nts/100) == 0:
        print 'timestep {0}/{1}, t={2}'.format(k, Nts, t)
        vfile = rdir + vfileprfx + '__t{0}.vtu'.format(t)
        pfile = rdir + pfileprfx + '__t{0}.vtu'.format(t)
        vu.writevp_paraview(velvec=old_v, pvec=p, vfile=vfile, pfile=pfile,
                            strtojson=visujsonstr)
        vfilelist.append(vfile)
        pfilelist.append(pfile)

vu.collect_vtu_files(vfilelist, vfileprfx+'.pvd')
vu.collect_vtu_files(pfilelist, pfileprfx+'.pvd')

print '\n### for visualization try:'
print 'paraview ' + vfileprfx + '.pvd'
print 'paraview ' + pfileprfx + '.pvd'

tikzfile = 'tikzs/p_drivcav-N{0}-tE{1}-Nts{2}-bccomg{3}'.\
    format(N, tE, Nts, omeg)
vu.plot_prs_outp(outsig=poutlist, tmesh=trange, fignum=222, tikzfile=tikzfile)

tikzfile = 'tikzs/v_drivcav-N{0}-tE{1}-Nts{2}-bccomg{3}'.\
    format(N, tE, Nts, omeg)
vu.plot_prs_outp(outsig=voutlist, tmesh=trange, fignum=123,
                 tikzfile=tikzfile)

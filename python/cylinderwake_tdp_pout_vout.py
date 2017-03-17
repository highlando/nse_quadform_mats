import numpy as np

import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla

import conv_tensor_utils as ctu
import visualization_utils as vu


N = 3
Re = 90

NVdict = {1: 5812, 3: 19468}

NV = NVdict[N]
savedmatsstr = '../data/cylinderwake__mats_NV{1}_Re{0}.mat'.\
    format(1, NV)
visujsonstr = '../data/visualization_cylinderwake_NV{0}.jsn'.format(NV)
print 'Re = ', Re
print 'NV = ', NV

# get the coefficients
mats = scipy.io.loadmat(savedmatsstr)
A = 1./Re*mats['A'] + mats['L1'] + mats['L2']
M = mats['M']
J = mats['J']
hmat = mats['H']
fv = mats['fv'] + 1./Re*mats['fv_diff'] + mats['fv_conv']
fp = mats['fp'] + mats['fp_div']
pcmat = mats['Cp']
vcmat = mats['Cv']
NV, NP = fv.shape[0], fp.shape[0]

t0, tE, Nts = 0., 4., 2**11
trange = np.linspace(t0, tE, Nts+1)
DT = (tE-t0)/Nts

rdir = 'results/'
vfileprfx = 'v_tdpcyl_NV{0}_Re{1}_tE{2}_Nts{3}'.format(NV, Re, tE, Nts)
pfileprfx = 'p_tdpcyl_NV{0}_Re{1}_tE{2}_Nts{3}'.format(NV, Re, tE, Nts)


sysmat = sps.vstack([sps.hstack([M+DT*A, -J.T]),
                     sps.hstack([J, sps.csc_matrix((NP, NP))])])

print 'computing Stokes solution to be used as initial value...'

fvstks = mats['fv'] + 1./Re*mats['fv_diff']
Astks = 1./Re*mats['A']

stksrhs = np.vstack([fvstks, fp])
stksmat = sps.vstack([sps.hstack([Astks, -J.T]),
                      sps.hstack([J, sps.csc_matrix((NP, NP))])])
stksvp = spsla.spsolve(stksmat, stksrhs).reshape((NV+NP, 1))
stksv = stksvp[:NV].reshape((NV, 1))
stksp = stksvp[NV:].reshape((NP, 1))

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
    crhsv = M*old_v + DT*(fv - ctu.eva_quadterm(hmat, old_v))
    crhs = np.vstack([crhsv, fp])
    vp_new = np.atleast_2d(sysmati(crhs.flatten())).T
    old_v = vp_new[:NV]
    p = vp_new[NV:]

    poutlist.append((pcmat*p)[0][0])
    voutlist.append((vcmat*old_v).flatten())
    if np.mod(k, Nts/10) == 0:
        print 'timestep {0}/{1}, t={2}, |v|={3}'.\
            format(k, Nts, t, np.linalg.norm(old_v))
        vfile = rdir + vfileprfx + '__t{0}.vtu'.format(t)
        pfile = rdir + pfileprfx + '__t{0}.vtu'.format(t)
        vu.writevp_paraview(velvec=old_v, pvec=p, vfile=vfile, pfile=pfile,
                            strtojson=visujsonstr)
        vfilelist.append(vfile)
        pfilelist.append(pfile)

vu.collect_vtu_files(vfilelist, vfileprfx+'.pvd')
vu.collect_vtu_files(pfilelist, pfileprfx+'.pvd')

tikzfile = 'tikzs/p_nsequadtens-N{0}-tE{1}-Nts{2}'.format(N, tE, Nts)
vu.plot_prs_outp(outsig=poutlist, tmesh=trange,
                 tikzfile=tikzfile)
tikzfile = 'tikzs/v_nsequadtens-N{0}-tE{1}-Nts{2}'.format(N, tE, Nts)
vu.plot_prs_outp(outsig=voutlist, tmesh=trange,
                 tikzfile=tikzfile)

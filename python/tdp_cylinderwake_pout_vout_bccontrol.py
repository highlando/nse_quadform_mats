import numpy as np

import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla

import conv_tensor_utils as ctu
import visualization_utils as vu

N = 2
Re = 40
palpha = 1e-3  # penalization parameter for the robin relaxation

NVdict = {2: 9384}

NV = NVdict[N]
savedmatsstr = '../data/cylinderwake_mats_' +\
    'N{1}_Re{0}_bccontrol_palpha{2}.mat'.format(1, NV, 1)

visujsonstr = '../data/visualization_cylN{0}.jsn'.format(N)

print 'Re = ', Re
print 'NV = ', NV

# get the coefficients
mats = scipy.io.loadmat(savedmatsstr)
A = 1./Re*mats['A'] + mats['L1'] + mats['L2'] + 1./palpha*mats['Arob']
Brob = 1./palpha*mats['Brob']
M = mats['M']
J = mats['J']
hmat = mats['H']
fv = mats['fv'] + 1./Re*mats['fv_diff'] + mats['fv_conv']
fp = mats['fp'] + mats['fp_div']
pcmat = mats['Cp']
vcmat = mats['Cv']
NV, NP = fv.shape[0], fp.shape[0]

t0, tE, Nts = 0., 6., 2**10
trange = np.linspace(t0, tE, Nts+1)
DT = (tE-t0)/Nts
omeg = 3.  # parameter for the frequency of the input signal

rdir = 'results/'
vfileprfx = 'v_tdpcyl_NV{0}_Re{1}_tE{2}_Nts{3}_bccomg{4}'.\
    format(NV, Re, tE, Nts, omeg)
pfileprfx = 'p_tdpcyl_NV{0}_Re{1}_tE{2}_Nts{3}_bccomg{4}'.\
    format(NV, Re, tE, Nts, omeg)


# define an input function
def bbcu(t):
    uvec = np.sin(omeg*np.pi*t/(tE-t0))*np.array([[1], [-1]])
    return np.dot(Brob, uvec)

sysmat = sps.vstack([sps.hstack([M+DT*A, -J.T]),
                     sps.hstack([J, sps.csc_matrix((NP, NP))])])

print 'computing Stokes solution to be used as initial value...'

fvstks = mats['fv'] + 1./Re*mats['fv_diff'] + bbcu(t0)
Astks = 1./Re*mats['A'] + 1./palpha*mats['Arob']

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
vu.cylwake_paraview(velvec=stksv, pvec=stksp, vfile=vfile, pfile=pfile,
                    strtojson=visujsonstr)
vfilelist, pfilelist = [vfile], [pfile]

for k, t in enumerate(trange):
    crhsv = M*old_v + DT*(fv - ctu.eva_quadterm(hmat, old_v) + bbcu(t))
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
        vu.cylwake_paraview(velvec=old_v, pvec=p, vfile=vfile, pfile=pfile,
                            strtojson=visujsonstr)
        vfilelist.append(vfile)
        pfilelist.append(pfile)

vu.collect_vtu_files(vfilelist, vfileprfx+'.pvd')
vu.collect_vtu_files(pfilelist, pfileprfx+'.pvd')

tikzfile = 'tikzs/p_nsequadtens-N{0}-tE{1}-Nts{2}-bccomg{3}'.\
    format(N, tE, Nts, omeg)
vu.plot_prs_outp(outsig=poutlist, tmesh=trange, tikzfile=tikzfile)

tikzfile = 'tikzs/v_nsequadtens-N{0}-tE{1}-Nts{2}-bccomg{3}'.\
    format(N, tE, Nts, omeg)
vu.plot_prs_outp(outsig=voutlist, tmesh=trange,
                 tikzfile=tikzfile)

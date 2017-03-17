import numpy as np
import scipy.io
import scipy.sparse as sps
import scipy.sparse.linalg as spsla
import conv_tensor_utils as ctu
import visualization_utils as vu


# hard coded paths and dictionary for data
NVdict          = {10: 722, 20: 3042,  30: 6962}
savedmatsstr    = lambda NV: '../data/drivencavity__mats_NV{1}_Re{0}.mat'.format(1,NV)
visujsonstr     = lambda N : '../data/visualization_drivencavity_N{0}.jsn'.format(N)


# setup parameters
Re          = 800
npicardstps = 5
N           = 20
NV          = NVdict[N]
omeg        = 4.  # parameter for the frequency of the input signal


# parameters for time stepping
t0          = 0.
tE          = 20.
Nts         = 2**12
DT          = (tE-t0)/Nts
trange      = np.linspace(t0, tE, Nts+1)


#parameters for results, directories
rdir        = 'results/'
vfileprfx   = 'v_tdpdrivcav_NV{0}_Re{1}_tE{2}_Nts{3}_bccomg{4}'.format(NV, Re, tE, Nts, omeg)
pfileprfx   = 'p_tdpdrivcav_NV{0}_Re{1}_tE{2}_Nts{3}_bccomg{4}'.format(NV, Re, tE, Nts, omeg)
poutlist    = []
voutlist    = []
vfile       = lambda t : rdir + vfileprfx + '__t{0}.vtu'.format(t)
pfile       = lambda t : rdir + pfileprfx + '__t{0}.vtu'.format(t)
vfilelist   = [vfile(trange[0])]
pfilelist   = [pfile(trange[0])]
ptikzfile   = 'tikz/p_nsequadtens-N{0}-tE{1}-Nts{2}-bccomg{3}'.format(N, tE, Nts, omeg)
vtikzfile   = 'tikz/v_nsequadtens-N{0}-tE{1}-Nts{2}-bccomg{3}'.format(N, tE, Nts, omeg)


# print reynolds number, discretization lvl, and other params
print 'Re  = {0}'.format(Re)
print 'NV  = {0}'.format(NV)
print 't0  = {0}'.format(t0)
print 'tE  = {0}'.format(tE)
print 'Nts = {0}'.format(Nts)
print 'DT  = {0}'.format(DT)
print '\n'


# load the coefficients matrices
mats    = scipy.io.loadmat(savedmatsstr(NV))
M       = mats['M']
A       = 1./Re*mats['A'] + mats['L1'] + mats['L2']
B       = mats['B']
J       = mats['J']
hmat    = mats['H']
fv      = mats['fv'] + 1./Re*mats['fv_diff'] + mats['fv_conv']
fp      = mats['fp'] + mats['fp_div']
pcmat   = mats['Cp']
vcmat   = mats['Cv']
NV, NP  = fv.shape[0], fp.shape[0]


# Fix the p
J  = J[:-1, :]
fp = fp[:-1, :]


# restrict it to two dofs in the input
NU = B.shape[1]
B  = B[:, [0, NU/2]]


# define an input function
def bbcu(t):
    uvec = np.array([[np.sin(omeg*np.pi*t/(tE-t0))],[np.cos(omeg*np.pi*t/(tE-t0))]])
    return B*uvec


# factorization of system matrix
print 'computing LU once...'
sysmat  = sps.vstack([sps.hstack([M+DT*A, -J.T]), sps.hstack([J, sps.csc_matrix((NP-1, NP-1))])]).tocsc()
sysmati = spsla.factorized(sysmat)


# compute stokes solution as initial value
print 'computing Stokes solution to be used as initial value...'
fvstks  = mats['fv'] + 1./Re*mats['fv_diff'] + bbcu(t0)
Astks   = 1./Re*mats['A']
stksrhs = np.vstack([fvstks, fp])
stksmat = sps.vstack([sps.hstack([Astks, -J.T]),sps.hstack([J, sps.csc_matrix((NP-1, NP-1))])]).tocsc()
stksvp  = spsla.spsolve(stksmat, stksrhs).reshape((NV+NP-1, 1))
stksv   = stksvp[:NV].reshape((NV, 1))
stksp   = (np.r_[stksvp[NV:].flatten(), 0]).reshape((NP, 1))



# Preparing for the output
vu.writevp_paraview(velvec=stksv, pvec=stksp, vfile=vfile(trange[0]), pfile=pfile(trange[0]), strtojson=visujsonstr(N))


# time stepping
print 'doing the time loop...'
old_v = stksv

for k, t in enumerate(trange):
    crhsv   = M*old_v + DT*(fv - ctu.eva_quadterm(hmat, old_v) + bbcu(t))
    crhs    = np.vstack([crhsv, fp])
    vp_new  = np.atleast_2d(sysmati(crhs.flatten())).T
    old_v   = vp_new[:NV]
    p       = (np.r_[vp_new[NV:].flatten(), 0]).reshape((NP, 1))

    poutlist.append((pcmat*p)[0][0])
    voutlist.append((vcmat*old_v).flatten())
    if np.mod(k, round(Nts/100)) == 0:
        print 'timestep {0:4d}/{1}, t={2:f}'.format(k, Nts, t)
        vu.writevp_paraview(velvec=old_v, pvec=p, vfile=vfile(t), pfile=pfile(t),strtojson=visujsonstr(N))
        vfilelist.append(vfile(t))
        pfilelist.append(pfile(t))


# save collection to pvd file
vu.collect_vtu_files(vfilelist, vfileprfx+'.pvd')
vu.collect_vtu_files(pfilelist, pfileprfx+'.pvd')


# write to tikz file
vu.plot_prs_outp(outsig=poutlist, tmesh=trange, fignum=222, tikzfile=ptikzfile)
vu.plot_prs_outp(outsig=voutlist, tmesh=trange, fignum=123, tikzfile=vtikzfile)



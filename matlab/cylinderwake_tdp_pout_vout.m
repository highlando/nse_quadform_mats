%% clear all
clearvars, close all, clc


%% hard coded paths and dictionary for data
NVdict          = [5812, 9356, 19468];
savedmatsstr    = @(NV) sprintf('%s/data/cylinderwake__mats_NV%d_Re%d.mat',fileparts(pwd),NV,1);
visujsonstr     = @(NV) sprintf('%s/data/visualization_cylinderwake_NV%d.jsn',fileparts(pwd),NV);


%% setup standard parameters
N           = 1;
Re          = 40;
npicardstps = 5;


%% parameters for time stepping
t0          = 0.;
tE          = 4.;
Nts         = 2^11;


%% further parameters
NV          = NVdict(N);
DT          = (tE-t0)/Nts;
trange      = linspace(t0, tE, Nts+1);


%% parameters for results, directories
rdir        = 'results/'
vfileprfx   = sprintf('v_tdpcyl_NV%d_Re%d_tE%e_Nts%d',NV, Re, tE, Nts);
pfileprfx   = sprintf('p_tdpcyl_NV%d_Re%d_tE%e_Nts%d'.NV, Re, tE, Nts);
poutlist    = {};
voutlist    = {};
vfile       = @(t) sprintf('%s%s__t%e.vtu',rdir,vfileprfx,t);
pfile       = @(t) sprintf('%s%s__t%e.vtu',rdir,pfileprfx,t);
vfilelist   = {vfile(trange(1))};
pfilelist   = {pfile(trange(1))};
ptikzfile   = sprintf('tikz/p_nsequadtens-N%d-tE%e-Nts%d', N, tE, Nts);
vtikzfile   = sprintf('tikz/v_nsequadtens-N%d-tE%e-Nts%d', N, tE, Nts);


%% print reynolds number and discretization lvl
fprintf('Re           = %d\n',Re);
fprintf('NV           = %d\n',NV);
fprintf('Picardsteps  = %d\n',npicardstps);
fprintf('t0           = %e\n',t0);
fprintf('tE           = %e\n',tE);
fprintf('DT           = %e\n',DT);
fprintf('\n')


%% load the coefficients matrices
mats    = load(savedmatsstr(NV));
M       = mats.M;
A       = 1./Re*mats.A + mats.L1 + mats.L2;
J       = mats.J;
hmat    = mats.H;
Brob    = mats.Brob;
fv      = mats.fv + 1./Re*mats.fv_diff + mats.fv_conv;
fp      = mats.fp + mats.fp_div;
NV      = size(fv,1);
NP      = size(fp,1);

# load the coefficients matrices
mats    = scipy.io.loadmat(savedmatsstr(NV))
M       = mats['M']
A       = 1./Re*mats['A'] + mats['L1'] + mats['L2']
J       = mats['J']
hmat    = mats['H']
fv      = mats['fv'] + 1./Re*mats['fv_diff'] + mats['fv_conv']
fp      = mats['fp'] + mats['fp_div']
pcmat   = mats['Cp']
vcmat   = mats['Cv']
NV, NP  = fv.shape[0], fp.shape[0]


# factorization of system matrix
print 'computing LU once...'
sysmat  = sps.vstack([sps.hstack([M+DT*A, -J.T]), sps.hstack([J, sps.csc_matrix((NP, NP))])]).tocsc()
sysmati = spsla.factorized(sysmat)


# compute stokes solution as initial value
print 'computing Stokes solution to be used as initial value...'
fvstks  = mats['fv'] + 1./Re*mats['fv_diff']
Astks   = 1./Re*mats['A']
stksrhs = np.vstack([fvstks, fp])
stksmat = sps.vstack([sps.hstack([Astks, -J.T]),sps.hstack([J, sps.csc_matrix((NP, NP))])]).tocsc()
stksvp  = spsla.spsolve(stksmat, stksrhs).reshape((NV+NP, 1))
stksv   = stksvp[:NV].reshape((NV, 1))
stksp   = stksvp[NV:].reshape((NP, 1))

# Preparing for the output
vu.writevp_paraview(velvec=stksv, pvec=stksp, vfile=vfile(trange[0]), pfile=pfile(trange[0]), strtojson=visujsonstr(NV))


# time stepping
print 'doing the time loop...'
old_v = stksv

for k, t in enumerate(trange):
    crhsv   = M*old_v + DT*(fv - ctu.eva_quadterm(hmat, old_v))
    crhs    = np.vstack([crhsv, fp])
    vp_new  = np.atleast_2d(sysmati(crhs.flatten())).T
    old_v   = vp_new[:NV]
    p       = vp_new[NV:]

    poutlist.append((pcmat*p)[0][0])
    voutlist.append((vcmat*old_v).flatten())
    if np.mod(k, round(Nts/10)) == 0:
        print('timestep {0:4d}/{1}, t={2:f}, |v|={3:e}'.format(k, Nts, t, np.linalg.norm(old_v)))
        vu.writevp_paraview(velvec=old_v, pvec=p, vfile=vfile(t), pfile=pfile(t),strtojson=visujsonstr(NV))
        vfilelist.append(vfile(t))
        pfilelist.append(pfile(t))


# save collection to pvd file
vu.collect_vtu_files(vfilelist, vfileprfx+'.pvd')
vu.collect_vtu_files(pfilelist, pfileprfx+'.pvd')

# write to tikz file
vu.plot_prs_outp(outsig=poutlist, tmesh=trange, tikzfile=ptikzfile)
vu.plot_prs_outp(outsig=voutlist, tmesh=trange, tikzfile=vtikzfile)





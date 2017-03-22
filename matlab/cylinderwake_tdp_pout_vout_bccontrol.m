function cylinderwake_tdp_pout_vout_bccontrol(N_, Re_, Picardsteps_, palpha_, omega_, t0_, tE_, Nts_)
%% CYLINDERWAKE_TDP_POUT_VOUT_BCCONTROL
%
%  Calling Sequences:
%
%  default  -   cylinderwake_tdp_pout_vout_bccontrol()
%               cylinderwake_tdp_pout_vout_bccontrol(N, Re, Picardsteps, palpha, omega, t0, tE, Nts)
%



%% hard coded paths and dictionary for data
NVdict          = [5824, 9384, 19512];
savedmatsstr    = @(NV) sprintf('%s/data/cylinderwake__mats_NV%d_Re%d_bccontrol_palpha%d.mat',fileparts(pwd),NV,1,1);
visujsonstr     = @(NV) sprintf('%s/data/visualization_cylinderwake_NV%d.jsn',fileparts(pwd),NV);


%% setup standard parameters
N           = 1;
Re          = 40;
npicardstps = 5;
palpha      = 1e-3;
omeg        = 3.;  % parameter for the frequency of the input signal


%% parameters for time stepping
t0          = 0.;
tE          = 4;
Nts         = 2^11;


%% get command line input and overwrite standard paramters if necessary
if nargin == 8
    N           = N_;
    Re          = Re_;
    npicardstps = Picardsteps_;
    palpha      = palpha_;
    omeg        = omega_;
    t0          = t0_;
    tE          = tE_;
    Nts         = Nts_;
elseif nargin ~=0
   error('Unkown Number of input arguments'); 
end


%% further parameters
NV          = NVdict(N);
DT          = (tE-t0)/Nts;
trange      = linspace(t0, tE, Nts+1);


%% parameters for results, directories
rdir        = 'results/';
vfileprfx   = sprintf('v_tdpcyl_NV%d_Re%d_tE%e_Nts%d_bccomg%e',NV, Re, tE, Nts, omeg);
pfileprfx   = sprintf('p_tdpcyl_NV%d_Re%d_tE%e_Nts%d_bccomg%e',NV, Re, tE, Nts, omeg);
voutlist    = cell(1,length(trange));
poutlist    = cell(1,length(trange));
vfile       = @(t) sprintf('%s%s__t%e.vtu',rdir,vfileprfx,t);
pfile       = @(t) sprintf('%s%s__t%e.vtu',rdir,pfileprfx,t);
vfilelist   = {vfile(trange(1))};
pfilelist   = {pfile(trange(1))};
vtikzfile   = sprintf('tikz/v_nsequadtens-N%d-tE%e-Nts%d_bccomg%e', N, tE, Nts, omeg);
ptikzfile   = sprintf('tikz/p_nsequadtens-N%d-tE%e-Nts%d_bccomg%e', N, tE, Nts, omeg);


%% print reynolds number and discretization lvl
fprintf('Re           = %d\n',Re);
fprintf('NV           = %d\n',NV);
fprintf('Picardsteps  = %d\n',npicardstps);
fprintf('palpha       = %e\n',palpha);
fprintf('omega        = %e\n',omeg);
fprintf('t0           = %e\n',t0);
fprintf('tE           = %e\n',tE);
fprintf('Nts          = %d\n',Nts);
fprintf('DT           = %e\n',DT);
fprintf('\n')


%% load the coefficients matrices
mats    = load(savedmatsstr(NV));
M       = mats.M;
A       = 1./Re*mats.A + mats.L1 + mats.L2 + 1./palpha*mats.Arob;
Brob    = 1./palpha*mats.Brob:
J       = mats.J;
hmat    = mats.H;
fv      = mats.fv + 1./Re*mats.fv_diff + mats.fv_conv;
fp      = mats.fp + mats.fp_div;
pcmat   = mats.Cp;
vcmat   = mats.Cv;
NV      = size(fv,1);
NP      = size(fp,1);


%% define an input function
bbcu = @(t) Brob*sin(omeg*pi*t/(tE-t0))*[1;-1];


%% factorization of system matrix
fprintf('computing LU once...\n');
sysmati                     = sparse([[M+DT*A, -J'];[J, sparse(NP, NP)]]);
[sysL,sysU,sysP,sysQ,sysR]  = lu(sysmati);            
sysmati                     = @(x) sysQ * (sysU \ (sysL \ (sysP * (sysR \ x))));


%% compute stokes solution as initial value
fprintf('computing Stokes solution to be used as initial value...\n');
fvstks  = mats.fv + 1./Re*mats.fv_diff + bbcu(t0);
Astks   = 1/Re*mats.A + 1/palpha*mats.Arob;
stksrhs = [fvstks; fp];
stksmat = sparse([[Astks, -J'];[J, sparse(NP, NP)]]);
stksvp  = stksmat\stksrhs;
stksv   = stksvp(1:NV);
stksp   = stksvp(NV+1:end);


%% Preparing for the output
writevp_paraview(stksv, stksp, visujsonstr(NV), vfile(trange(1)), pfile(trange(1)));


%% time stepping
fprintf('doing the time loop...\n');
old_v = stksv;

for k = 1:length(trange)
    crhsv   = M*old_v + DT*(fv - eva_quadterm(hmat, old_v) + bbcu(t));
    crhs    = [crhsv; fp];
    vp_new  = sysmati(crhs);
    old_v   = vp_new(1:NV);
    p       = vp_new(NV+1:end);
    
    poutlist{ k } = pcmat*p;
    voutlist{ k } = vcmat*old_v;
    if mod(k, round(Nts/10)) == 0
        fprintf('timestep %4d/%d, t=%f, |v|=%e\n',k, Nts, trange(k), norm(old_v));
        writevp_paraview(old_v, p, visujsonstr(NV), vfile(trange(k)), pfile(trange(k)));
        pfilelist{ length(pfilelist) + 1 } = pfile(trange(k));
        vfilelist{ length(vfilelist) + 1 } = vfile(trange(k));
    end
    
end


%% save collection to pvd file
collect_vtu_files(vfilelist, strcat(vfileprfx,'.pvd'));
collect_vtu_files(pfilelist, strcat(pfileprfx,'.pvd'));


%vu.plot_prs_outp(outsig=poutlist, tmesh=trange, tikzfile=ptikzfile)
%vu.plot_prs_outp(outsig=voutlist, tmesh=trange, tikzfile=vtikzfile)


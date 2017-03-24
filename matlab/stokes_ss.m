%% clear all
clear all, close all, clc


%% hard coded paths and dictionary for data
NVdict          = [5812, 9356, 19468];
savedmatsstr    = @(NV) sprintf('%s/data/cylinderwake__mats_NV%d_Re%d.mat',fileparts(pwd),NV,1);
visujsonstr     = @(NV) sprintf('%s/data/visualization_cylinderwake_NV%d.jsn',fileparts(pwd),NV);


%% setup standard parameters
N           = 1;
Re          = 80;
npicardstps = 5;


%% visualisation files
NV    = NVdict(N);
pfile = sprintf('p__cylinderwake_stst_Re%d_NV%d.vtu',Re, NV);
vfile = sprintf('v__cylinderwake_stst_Re%d_NV%d.vtu',Re, NV);


%% print reynolds number and discretization lvl
fprintf('Re           = %d\n',Re);
fprintf('NV           = %d\n',NV);
fprintf('Picardsteps  = %d\n',npicardstps);
fprintf('pfile        = %s\n',pfile);
fprintf('vfile        = %s\n',vfile);
fprintf('\n')


%% load the coefficients matrices
mats    = load(savedmatsstr(NV));
M       = mats.M;
A       = 1./Re*mats.A;
J       = mats.J;
hmat    = mats.H;
fv      = mats.fv + 1./Re*mats.fv_diff;
fp      = mats.fp + mats.fp_div;
NV      = size(fv,1);
NP      = size(fp,1);


%% solve stokes equation
stksrhs   = [fv; fp];
stksmat   = sparse([[A, -J'];[J, sparse(NP, NP)]]);
stksvp    = stksmat\stksrhs;
stksv     = stksvp(1:NV);
stksp     = stksvp(NV+1:end);


%% print results
resstksmom  = A*stksv - J'*stksp - fv;
resconti    = J*stksv - fp;

fprintf('*** Done ***\n');
fprintf('The Stokes momentum eq residual: %e\n',norm(resstksmom));
fprintf('The conti residual:              %e\n\n',norm(resconti));

vssstks = mats.v_ss_stokes;
pssstks = mats.p_ss_stokes;
 
fprintf('Difference of v solutions:       %e\n',(norm(stksv - vssstks)));
fprintf('Difference of p solutions:       %e\n',(norm(stksp - pssstks)));



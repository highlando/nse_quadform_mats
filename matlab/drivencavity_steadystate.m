function drivencavity_steadystate(N_,Re_,Picardsteps_)
%% DRIVENCAVITY_STEADYSTATE
%
%  Calling Sequences:
%
%  default  -   drivencavity_steadystate()
%               drivencavity_steadystate(N, Re, Picardsteps)
%


%% hard coded paths and dictionary for data
NVdict          = [722, 3042, 6962];
savedmatsstr    = @(NV) sprintf('%s/data/drivencavity__mats_NV%d_Re%d.mat',fileparts(pwd),NV,1);
visujsonstr     = @(NV) sprintf('%s/data/visualization_drivencavity_NV%d.jsn',fileparts(pwd),NV);


%% setup standard parameters
N           = 1;
Re          = 40;
npicardstps = 5;


%% get command line input and overwrite standard paramters if necessary
if nargin == 3
    N           = N_;
    Re          = Re_;
    npicardstps = Picardsteps_;
elseif nargin ~=0
   error('Unkown Number of input arguments'); 
end

%% visualisation files
NV    = NVdict(N);
pfile = sprintf('results/p__drivencavity_stst_Re%d_NV%d.vtu',Re, NV);
vfile = sprintf('results/v__drivencavity_stst_Re%d_NV%d.vtu',Re, NV);


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
A       = 1./Re*mats.A + mats.L1 + mats.L2;
J       = mats.J;
hmat    = mats.H;
fv      = mats.fv + 1./Re*mats.fv_diff + mats.fv_conv;
fp      = mats.fp + mats.fp_div;
NV      = size(fv,1);
NP      = size(fp,1);


%% solve steady state equation
updnorm     = 1;
curv        = zeros(NV, 1);
curvp       = zeros(NV+NP, 1);
stpcount    = 0;

while updnorm > 1e-10
    
    picard = (stpcount < npicardstps);
    
    [H1k, H2k] = linearzd_quadterm(hmat, curv);
    if picard
        currhs = [fv; fp];
        HL     = H1k;
    else
        currhs = [fv+eva_quadterm(hmat, curv); fp];
        HL     = H1k+H2k;
    end
    
    cursysmat = sparse([[A+HL, -J'];[J, sparse(NP, NP)]]);
    nextvp    = cursysmat\currhs;
    
    if picard
        fprintf('Iteration step %d (Picard)\n',stpcount);
    else
        fprintf('Iteration step %d (Newton)\n',stpcount);
    end
    
    nextv       = nextvp(1:NV);
    nextp       = nextvp(NV+1:end);
    
    curnseres   = A*nextv + eva_quadterm(hmat,nextv) - J'*nextp - fv;
    fprintf('Norm of nse residual:   %e\n',norm(curnseres));
    updnorm = norm(nextv - curv) / norm(nextv);
    fprintf('Norm of current update: %e\n\n\n',norm(updnorm));
    curv = nextv;
    stpcount = stpcount + 1;
end


%% print results
fprintf('*** Done ***\n');
fprintf('NSE momentum eq residual: %e\n',norm(curnseres));
resconti = J*nextv - fp;
fprintf('The conti residual:       %e\n\n',norm(resconti));


%% write paraview
writevp_paraview(nextv, nextp, visujsonstr(NV), vfile, pfile);
fprintf('*** for visualization try ***\n');
fprintf('paraview %s\n',vfile);
fprintf('paraview %s\n',pfile);

function [L1,L2] = linearzd_quadterm(H, linv)
%% LINEARZED_QUADTERM - compute the matrices L1, L2 that represent the linearized convection
%
%    H(v, v) ~ L1*v + L2*v - H(linv, linv)
% 
% 
%     Parameters:
%     
%     H : (nv, nv*nv) sparse array
%         the tensor (as a matrix) that evaluates the convection term
%     linv : (nv, 1)  array
%         the stat at which the linearization is about 
% 
%     


%% check out arguments
if nargout < 1 || nargout > 2 
   error('Wrong number of output arguments.\n');
end

%% perform operation
fprintf('assembling hlmat ...\n');
nv = size(linv,1);

if nargout == 2
    L1 = H * sparse((kron(speye(nv), linv)));
    L2 = H * sparse((kron(linv, speye(nv))));
else
    L1 = H * sparse((kron(speye(nv), linv) + kron(linv, speye(nv))));
end
 




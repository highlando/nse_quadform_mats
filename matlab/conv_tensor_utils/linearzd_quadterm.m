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


%% check input arguments
assert(isvector(linv),'linv is no vector.');


%% perform operation
fprintf('assembling hlmat ...\n');
nv = size(linv,1);


%% perform L1 = H * sparse((kron(speye(nv), linv)));
slinv = sparse(linv);
for k = 1:nv
   L1(:,k)=H(:,(k-1)*nv+1:k*nv)*slinv; 
end


%% perform L2 = H * sparse((kron(linv, speye(nv))));
L2 = sparse(nv,nv);
for k = 1:nv
   L2 = L2 + linv(k)*H(:,(k-1)*nv+1:k*nv);
end


%% return results
if nargout == 1
    L1 = L1 + L2;
end
 




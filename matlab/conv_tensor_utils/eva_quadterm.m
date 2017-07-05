function hvv =  eva_quadterm(H, v)
%% EVA_QUADTERM function to evaluate H*kron(v, v) without forming kron(v, v)
%
%    Parameters:
%
%    H : (nv, nv*nv) sparse array
%        the tensor (as a matrix) that evaluates the convection term
%    v : (nv,1) matrix for the product
%
%

NV = size(v,1);
hvv = zeros(NV, 1);

for k = 1:length(v)
    hvv = hvv + H(:,(k-1)*NV+1:k*NV)*(v(k)*v);
end
   
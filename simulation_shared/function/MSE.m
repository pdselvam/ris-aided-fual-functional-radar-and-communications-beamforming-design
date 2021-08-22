function [e] = MSE(para, c, P, g)
%The Mean Square Error at the MMSE receiber
%  [e] = MSE(para, c, r, P, Rq)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   P: beamformer at BS
%   g: MMSE receiver
%Outputs:
%   e: mean square error
%Date: 14/07/2021
%Author: Zhaolin Wang

% MSE at the MMSE receiver
e = zeros(para.K,1);
for k = 1:para.K
    pk = P(:,k);
    ck = c(:,k);
    gk = g(k);
    
    e(k) = pow_abs(gk, 2) * ( sum(pow_abs(ck'*P,2)) + para.n )...
        - 2*real( gk*ck'*pk ) + 1;
end


end


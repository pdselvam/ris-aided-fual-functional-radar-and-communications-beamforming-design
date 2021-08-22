function [e] = MSE(para, c, r, P, Rq, g)
%The Mean Square Error at the MMSE receiber
%  [e] = MSE(para, c, r, P, Rq)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   r: equivalent radar channel
%   P: beamformer at BS
%   Rq: covariance matrix of radar signal
%   g: MMSE receiver
%Outputs:
%   e: mean square error
%Date: 30/05/2021
%Author: Zhaolin Wang

% MSE at the MMSE receiver
e = [];
for k = 1:para.K
    pk = P(:,k);
    ck = c(:,k);
    rk = r(:,k);
    gk = g(k);
    
    ek = pow_abs(gk, 2) * ( sum(pow_abs(ck'*P,2)) + real(rk'*Rq*rk) + para.n )...
        - 2*real( gk*ck'*pk ) + 1;
    e = [e;ek];
end


end


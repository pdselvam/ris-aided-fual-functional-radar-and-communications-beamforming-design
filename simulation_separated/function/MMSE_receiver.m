function [g] = MMSE_receiver(para, c, r, P, Rq)
%The MMSE receiver
%  [g] = MMSE_receiver(para, c, r, P, Rq)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   r: equivalent radar channel
%   P: beamformer at BS
%   Rq: covariance matrix of radar signal
%Outputs:
%   g: MMSE receiver
%Date: 30/05/2021
%Author: Zhaolin Wang

g = zeros(para.K,1);
for k = 1:para.K
    pk = P(:,k);
    ck = c(:,k);
    rk = r(:,k);
    g(k) = pk'*ck / ( sum(pow_abs(ck'*P,2)) + real(rk'*Rq*rk) + para.n );  
end

end


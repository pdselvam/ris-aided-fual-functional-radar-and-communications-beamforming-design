function [gamma] = SINR(para, c, r, P, Rq)
%Calculate the Signal to Interference and Noise Ratio
%  [gamma] = SINR(para, c, r, P, Rq)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   r: equivalent radar channel
%   P: beamformer at BS
%   Rq: covariance matrix of radar signal
%Outputs:
%   gamma: SINR
%Date: 30/05/2021
%Author: Zhaolin Wang

gamma = zeros(para.K,1);
for k = 1:para.K
    pk = P(:,k);
    P_inter = P; P_inter(:,k) = [];
    ck = c(:,k);
    rk = r(:,k);
    
    % SINR
    gamma(k) = pow_abs(ck'*pk,2) / ( para.n + sum(pow_abs(ck'*P_inter,2)) + real(rk'*Rq*rk) );
end
end


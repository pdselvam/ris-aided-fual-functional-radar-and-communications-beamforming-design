function [gamma] = SINR(para, c, P)
%Calculate the Signal to Interference and Noise Ratio
%  [gamma] = SINR(para, c, P)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   P: beamformer at BS
%Outputs:
%   gamma: SINR
%Date: 14/07/2021
%Author: Zhaolin Wang

gamma = zeros(para.K,1);
for k = 1:para.K
    pk = P(:,k);
    P_inter = P; P_inter(:,k) = [];
    ck = c(:,k);
    
    % SINR
    gamma(k) = pow_abs(ck'*pk,2) / ( para.n + sum(pow_abs(ck'*P_inter,2)));
end
end


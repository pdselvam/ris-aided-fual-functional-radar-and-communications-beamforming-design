function [g] = MMSE_receiver(para, c, P)
%The MMSE receiver
%  [g] = MMSE_receiver(para, c, r, P, Rq)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   P: beamformer at BS
%Outputs:
%   g: MMSE receiver
%Date: 14/07/2021
%Author: Zhaolin Wang

g = zeros(para.K,1);
for k = 1:para.K
    pk = P(:,k);
    ck = c(:,k);
    g(k) = pk'*ck / ( sum(pow_abs(ck'*P,2)) + para.n );  
end

end


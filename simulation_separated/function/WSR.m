function [wsr] = WSR(para, c, r, P, Rq)
%The Weighted Sum Rate of the system
%  [wsr] = WSR(para, c, r, P, Rq)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   r: equivalent radar channel
%   P: beamformer at BS
%   Rq: covariance matrix of radar signal
%Outputs:
%   wsr: weighted sum rate
%Date: 30/05/2021
%Author: Zhaolin Wang

gamma = SINR(para, c, r, P, Rq);
R = log2(1+gamma);

wsr = para.weight' * R;

end


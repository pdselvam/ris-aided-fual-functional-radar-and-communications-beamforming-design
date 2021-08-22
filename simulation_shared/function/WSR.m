function [wsr] = WSR(para, c, P)
%The Weighted Sum Rate of the system
%  [wsr] = WSR(para, c, P)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   r: equivalent radar channel
%   P: beamformer at BS
%Outputs:
%   wsr: weighted sum rate
%Date: 14/07/2021
%Author: Zhaolin Wang

gamma = SINR(para, c, P);
R = log2(1+gamma);

wsr = para.weight' * R;

end


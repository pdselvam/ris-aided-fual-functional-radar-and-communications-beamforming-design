function [pattern] = beampattern(P,Rq,Mc,Mr,theta_degree)
%The beampattern of ULA in a DFRC BS with seperate depolyment
%  [pattern] = beampattern(P,Rx,Mc,Mr,theta_degree)
%Inputs:
%   P: beamforming for communication signal
%   Rq: covariance matrix of the radar signal
%   Mc: number of communication antennas
%   Mr: number of radar antennas
%   theta_degree: range of direction angle
%Outputs:
%   pattern: beampattern
%Date: 28/02/2021
%Author: Zhaolin Wang

theta = theta_degree*pi/180;
pattern = zeros(length(theta),1);
% C = [Rx, zeros(Ntr, Ntc); zeros(Ntr,Ntc), P*P'];

for i = 1:length(theta)
    t = theta(i);
    a = ULA_func(t,Mc+Mr);
    a1 = a(1:Mr);
    a2 = a(Mr+1:end);
    pattern(i) = real(a1'*Rq*a1 + a2'*(P*P')*a2);
end

pattern = 10*log10(pattern);
end


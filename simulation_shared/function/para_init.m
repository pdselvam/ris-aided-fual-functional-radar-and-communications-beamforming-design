function [values] = para_init()
%Construct a struct of the initial values for all the parameters 
%  [values] = para_init()
%Inputs:
%   None
%Outputs:
%   values: a struct
%Date: 28/05/2021
%Author: Zhaolin Wang

values.noise = -117; % noise powe in dBm

values.M = 16; % overall antennas
% values.N = 20; % reflecting elements at RIS
values.N = 20; % reflecting elements at RIS

values.Pt = 10^(20/10); % overall transmit power
values.n = 1; % equivalent noise power
values.K = 4; % user number
values.phi_m = 0; % desired direction
values.weight = ones(values.K,1) ./ values.K; % weight of weighted sum rate 
values.pathloss_indirect = @(d) 35.6 + 22*log10(d); % path loss with d in m
values.pathloss_direct =  @(d) 32.6 + 36.7*log10(d); % path loss with d in m

values.BS_loc = [0,0]; % location of BS
values.RIS_loc = [200,0]; % location of RIS in m
values.user_center = [200, 30];
values.user_range = [0, 10]; % range of user location from RIS in m

values.rician = 1000; % rician factor
values.G = 1; % number of groups in IRS model

end





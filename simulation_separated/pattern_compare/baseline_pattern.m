clc
clear all
close all

rho = 500;


addpath('./function/')
%% Parameters
para = para_init();
d_BR = sqrt(para.RIS_loc(1)^2 + para.RIS_loc(2)^2);
theta_degree = -90:90;

%% Steering Vector
a = ULA_func(para.phi_m,para.M);
ar = a(1:para.Mr);
ac = a(para.Mr+1:end);
Z = para.Mc*eye(para.Mc) - ac*ac';

%% Generate user location
[user_loc, angle.user, d_RU, d_BU] = generate_user_location(para);
angle.RIS = atan(para.RIS_loc(2)/para.RIS_loc(1)); % direction of RIS from BS
angle.BS = pi + angle.RIS - pi/4; % direction of BS from RIS
% load('user_location.mat');

%% Path loss
path_loss.BU = para.pathloss_direct(d_BU)';
path_loss.BRU = para.pathloss_indirect(d_BR) + para.pathloss_indirect(d_RU)';
path_loss.BU = sqrt(10.^((-para.noise-path_loss.BU)/10));
path_loss.BRU = sqrt(10.^((-para.noise-path_loss.BRU)/10));


%% Monte Carlo Simulation
ite = 1;
wsr_all = zeros(ite,1);
prob_power_all = zeros(ite,1);
pattern_all = zeros(length(theta_degree),ite);

for step = 1:ite
    % Channel
    load('channel_pattern.mat');
    
    % Optimization Algorithm
    [P, Rq, wsr] = WMMSE_SDP(para, dc, dr, ar, Z, rho);
    
    % Beampattern
    [pattern] = beampattern(P, Rq, para.Mc, para.Mr, theta_degree);
    
    % Probing power
    prob_power = real(ar'*Rq*ar + ac'*(P*P')*ac);
    prob_power = 10*log10(prob_power);
    
    
    wsr_all(step) = wsr;
    prob_power_all(step) = prob_power;
    pattern_all(:,step) = pattern;
end

wsr_average = mean(wsr_all);
prob_power_average = mean(prob_power_all);
pattern_average = mean(pattern_all,2);

save('baseline_data_pattern.mat','para','wsr_average','prob_power_average','pattern_average');

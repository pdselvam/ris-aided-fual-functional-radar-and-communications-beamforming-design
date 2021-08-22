clc
clear all
close all

rho = 0.00000001;

addpath('../function/')
%% Parameters
para = para_init();
d_BR = sqrt(para.RIS_loc(1)^2 + para.RIS_loc(2)^2);
theta_degree = -90:90;

%% Steering Vector
a = ULA_func(para.phi_m,para.M);
Z = para.M*eye(para.M) - a*a';

%% Generate user location
[user_loc, angle.user, d_RU, d_BU] = generate_user_location(para);
angle.RIS = atan(para.RIS_loc(2)/para.RIS_loc(1)); % direction of RIS from BS
angle.BS = pi + angle.RIS - pi/4; % direction of BS from RIS
% plot_location(para, user_loc);
% save('user_location.mat','user_loc','angle','d_RU','d_BU');
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
%     [h, H, d, ~, ~] = generate_channel(para, angle, path_loss, a);
    load('channel_pattern.mat');
    
    % Optimization Algorithm
    [P, Theta, wsr] = WMMSE_FP(para, h, H, d, Z, rho);
    
    % Beampattern
    [pattern] = beampattern(P, para.M, theta_degree);
    
    % Probing power
    prob_power = real(a'*(P*P')*a);
    prob_power = 10*log10(prob_power);
    
    wsr_all(step) = wsr;
    prob_power_all(step) = prob_power;
    pattern_all(:,step) = pattern;
end

wsr_average = mean(wsr_all);
prob_power_average = mean(prob_power_all);
pattern_average = mean(pattern_all,2);

save('ris_data_single_pattern.mat','para','wsr_average','prob_power_average','pattern_average');

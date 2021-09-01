clc
clear all
close all


addpath('./function/')
%% Parameters
para = para_init();
d_BR = sqrt(para.RIS_loc(1)^2 + para.RIS_loc(2)^2);
theta_degree = -90:90;
rho_all = [10,100,1000,100000];

%% Steering Vector
a = ULA_func(para.phi_m,para.M);
Z = para.M*eye(para.M) - a*a';

%% Generate user location
load('user_location.mat');

%% Path loss
path_loss.BU = para.pathloss_direct(d_BU)';
path_loss.BRU = para.pathloss_indirect(d_BR) + para.pathloss_indirect(d_RU)';
path_loss.BU = sqrt(10.^((-para.noise-path_loss.BU)/10));
path_loss.BRU = sqrt(10.^((-para.noise-path_loss.BRU)/10));


% Channel
load('channel.mat');
figure();
for i = 1:length(rho_all)
    rho = rho_all(i);
    % Optimization Algorithm
    [~, ~, ~, wsr_all] = WMMSE_FP(para, h, H, d, Z, rho);

    len = length(wsr_all);
    if len < 30
        wsr_all = [wsr_all,wsr_all(end)*ones(1, 30-len)];
    else
        wsr_all = wsr_all(1:30);
    end
    plot(1:30,wsr_all,'LineWidth',1);
    hold on;
end

grid on;
legend('$\rho=10$','$\rho=100$','$\rho=1000$','$\rho=100000$','FontSize',12,'interpreter','latex');
ylabel('WSR [bps/Hz]','FontSize',12,'interpreter','latex');
xlabel('Steps','FontSize',12,'interpreter','latex');
xlim([1,30])
    
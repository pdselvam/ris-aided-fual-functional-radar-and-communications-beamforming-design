clc
clear all
close all

addpath('./function/')
%% trade-off
figure;
load('baseline.mat');
wsr_baseline = [0; wsr_baseline;wsr_baseline(end)];
prob_power_baseline = [prob_power_baseline(1); prob_power_baseline;0];
plot(wsr_baseline, prob_power_baseline,'-+b','LineWidth',3);
hold on;

load('./ris_aided_data/rician_0/ris_aided_single.mat');
pattern_ris_single = pattern_ris;
wsr_ris = [0; wsr_ris; wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-+g','LineWidth',3);

load('./ris_aided_data/rician_0/ris_aided_fully.mat');
pattern_ris_fully = pattern_ris;
wsr_ris = [0; wsr_ris;wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-+r','LineWidth',3);

grid on;
legend('without RIS','RIS-aided (single)','RIS-aided (fully)','FontSize',12,'interpreter','latex');
xlabel('WSR [bps/Hz]','FontSize',12,'interpreter','latex');
ylabel('Probing Power [dBm]','FontSize',12,'interpreter','latex');
ylim([19 34]);
xlim([0 10]);

%% radar only beampattern
para = para_init();
theta_degree = -90:90;
theta = theta_degree*pi/180;

pattern_radar_only = zeros(length(theta),1);
a = ULA_func(0,para.M);
R = para.Pt/para.M*(a*a');

for i = 1:length(theta)
    t = theta(i);
    a = ULA_func(t,para.M);
    pattern_radar_only(i) = real(a'*R*a);
end

pattern_radar_only = 10*log10(pattern_radar_only);

%% beampattern WSR = 7.9
figure;
plot(theta_degree, pattern_radar_only,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);
hold on;

% load('./pattern_data/rician_0/high_rate/ris_data_pattern.mat');
% plot(theta_degree, pattern_average,'-r','LineWidth',2);
% load('./pattern_data/rician_0/high_rate/ris_data_single_pattern.mat');
% plot(theta_degree, pattern_average,'-.g','LineWidth',2);
% load('./pattern_data/rician_0/high_rate/baseline_data_pattern.mat');
% plot(theta_degree, pattern_average,':b','LineWidth',2);
plot(theta_degree, pattern_ris_fully(:,11),'-r','LineWidth',2);
plot(theta_degree, pattern_ris_single(:,13),'-.g','LineWidth',2);
plot(theta_degree, pattern_baseline(:,16),':b','LineWidth',2);
grid on;

ylim([0,35]);
xlim([-90,90]);

legend('Radar-only','RIS-aided (fully)','RIS-aided (single)','without RIS','FontSize',12,'interpreter','latex');
xlabel('Azimuth Angle [degree]','FontSize',12,'interpreter','latex');
ylabel('Beampattern Value [dBm]','FontSize',12,'interpreter','latex');


%% beampattern WSR = 1.8
figure;
plot(theta_degree, pattern_radar_only,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);
hold on;
% 
% load('./pattern_data/rician_0/low_rate/ris_data_pattern.mat');
% plot(theta_degree, pattern_average,'-r','LineWidth',2);
% hold on;
% load('./pattern_data/rician_0/low_rate/ris_data_single_pattern.mat');
% plot(theta_degree, pattern_average,'-.g','LineWidth',2);
% load('./pattern_data/rician_0/low_rate/baseline_data_pattern.mat');
% plot(theta_degree, pattern_average,':b','LineWidth',2);
plot(theta_degree, pattern_ris_fully(:,2),'-r','LineWidth',2);
plot(theta_degree, pattern_ris_single(:,1),'-.g','LineWidth',2);
plot(theta_degree, pattern_baseline(:,3),':b','LineWidth',2);
grid on;

ylim([0,35]);
xlim([-90,90]);

legend('Radar-only','RIS-aided (fully)','RIS-aided (single)','without RIS','FontSize',12,'interpreter','latex');
xlabel('Azimuth Angle [degree]','FontSize',12,'interpreter','latex');
ylabel('Beampattern Value [dBm]','FontSize',12,'interpreter','latex');
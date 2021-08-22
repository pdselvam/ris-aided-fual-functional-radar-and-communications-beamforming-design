clc
clear all
close all
addpath('./simulation_separated/function/');
%% trade-off
%% baseline
figure;
load('./simulation_separated/baseline.mat');
pattern_baseline_sep = pattern_baseline;
wsr_baseline = [0; wsr_baseline;wsr_baseline(end)];
prob_power_baseline = [prob_power_baseline(1); prob_power_baseline;0];
plot(wsr_baseline, prob_power_baseline,'-^b','LineWidth',3,'MarkerSize',3);
hold on;

load('./simulation_shared/baseline.mat');
pattern_baseline_sh = pattern_baseline;
wsr_baseline = [0; wsr_baseline;wsr_baseline(end)];
prob_power_baseline = [prob_power_baseline(1); prob_power_baseline;0];
plot(wsr_baseline, prob_power_baseline,'-+b','LineWidth',3);

%% 20 elements
load('./simulation_separated/ris_aided_data/rician_1000/ris_aided_single.mat');
pattern_ris_20_sep = pattern_ris;
wsr_ris = [0; wsr_ris; wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-^g','LineWidth',3,'MarkerSize',3);

load('./simulation_shared/ris_aided_data/rician_1000/ris_aided_single.mat');
pattern_ris_20_sh = pattern_ris;
wsr_ris = [0; wsr_ris;wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-+g','LineWidth',3);

%% 60 elements
load('./simulation_separated/ris_aided_data/rician_1000/ris_aided_single_60.mat');
pattern_ris_60_sep = pattern_ris;
wsr_ris = [0; wsr_ris; wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-^c','LineWidth',3,'MarkerSize',3);

load('./simulation_shared/ris_aided_data/rician_1000/ris_aided_single_60.mat');
pattern_ris_60_sh = pattern_ris;
wsr_ris = [0; wsr_ris;wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-+c','LineWidth',3);

%% 100 elements
load('./simulation_separated/ris_aided_data/rician_1000/ris_aided_single_100.mat');
pattern_ris_100_sep = pattern_ris;
wsr_ris = [0; wsr_ris;wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-^r','LineWidth',3,'MarkerSize',3);

load('./simulation_shared/ris_aided_data/rician_1000/ris_aided_single_100.mat');
pattern_ris_100_sh = pattern_ris;
wsr_ris = [0; wsr_ris; wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-+r','LineWidth',3);

grid on;
legend('without RIS(separated)','without RIS(shared)','20 elements(separeted)','20 elements(shared)','60 elements(separated)','60 elements(shared)','100 elements(separated)','100 elements(shared)','FontSize',12,'interpreter','latex');
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


%% beampattern WSR = 5.4;
%% baseline
figure;
plot(theta_degree, pattern_radar_only,':','LineWidth',2,'Color',[0.5 0.5 0.5]);
hold on;
plot(theta_degree, pattern_baseline_sep(:,16),'-b','LineWidth',2);
plot(theta_degree, pattern_baseline_sh(:,8),'-r','LineWidth',2);
ylim([0,35]);
xlim([-90,90]);
grid on;
legend('Radar-only','without RIS(separated)','without RIS(shared)','FontSize',12,'interpreter','latex');
xlabel('Azimuth Angle [degree]','FontSize',12,'interpreter','latex');
ylabel('Beampattern Value [dBm]','FontSize',12,'interpreter','latex');

%% 20 elements
figure;
plot(theta_degree, pattern_radar_only,':','LineWidth',2,'Color',[0.5 0.5 0.5]);
hold on;
plot(theta_degree, pattern_ris_20_sep(:,14),'-b','LineWidth',2);
plot(theta_degree, pattern_ris_20_sh(:,8),'-r','LineWidth',2);
ylim([0,35]);
xlim([-90,90]);
grid on;
legend('Radar-only','20 elements(separated)','20 elements(shared)','FontSize',12,'interpreter','latex');
xlabel('Azimuth Angle [degree]','FontSize',12,'interpreter','latex');
ylabel('Beampattern Value [dBm]','FontSize',12,'interpreter','latex');

%% 60 elements
figure;
plot(theta_degree, pattern_radar_only,':','LineWidth',2,'Color',[0.5 0.5 0.5]);
hold on;
plot(theta_degree, pattern_ris_60_sep(:,12),'-b','LineWidth',2);
plot(theta_degree, pattern_ris_60_sh(:,5),'-r','LineWidth',2);
ylim([0,35]);
xlim([-90,90]);
grid on;
legend('Radar-only','60 elements(separated)','60 elements(shared)','FontSize',12,'interpreter','latex');
xlabel('Azimuth Angle [degree]','FontSize',12,'interpreter','latex');
ylabel('Beampattern Value [dBm]','FontSize',12,'interpreter','latex');

%% 100 elements
figure;
plot(theta_degree, pattern_radar_only,':','LineWidth',2,'Color',[0.5 0.5 0.5]);
hold on;
plot(theta_degree, pattern_ris_100_sep(:,11),'-b','LineWidth',2);
plot(theta_degree, pattern_ris_100_sh(:,4),'-r','LineWidth',2);
ylim([0,35]);
xlim([-90,90]);
grid on;
legend('Radar-only','100 elements(separated)','100 elements(shared)','FontSize',12,'interpreter','latex');
xlabel('Azimuth Angle [degree]','FontSize',12,'interpreter','latex');
ylabel('Beampattern Value [dBm]','FontSize',12,'interpreter','latex');

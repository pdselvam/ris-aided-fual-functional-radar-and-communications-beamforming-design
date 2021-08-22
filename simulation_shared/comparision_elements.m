clc
clear all
close all
addpath('./function/');

%% trade-off
load('baseline.mat');
wsr_baseline = [0; wsr_baseline;wsr_baseline(end)];
prob_power_baseline = [prob_power_baseline(1); prob_power_baseline;0];
plot(wsr_baseline, prob_power_baseline,'-+b','LineWidth',3);
hold on;


load('./ris_aided_data/rician_1000/ris_aided_single.mat');
pattern_ris_20 = pattern_ris;
wsr_ris = [0; wsr_ris;wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-+g','LineWidth',3);
hold on;


load('./ris_aided_data/rician_1000/ris_aided_single_60.mat');
pattern_ris_60 = pattern_ris;
wsr_ris = [0; wsr_ris;wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-+c','LineWidth',3);

load('./ris_aided_data/rician_1000/ris_aided_single_100.mat');
pattern_ris_100 = pattern_ris;
wsr_ris = [0; wsr_ris; wsr_ris(end)];
prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
plot(wsr_ris, prob_power_ris,'-+r','LineWidth',3);

grid on;
legend('without RIS','20 elements','60 elements','100 elements','FontSize',12,'interpreter','latex');
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
plot(theta_degree, pattern_baseline(:,16),':b','LineWidth',2);
plot(theta_degree, pattern_ris_20(:,13),'-.g','LineWidth',2);
plot(theta_degree, pattern_ris_60(:,10),'--c','LineWidth',2);
plot(theta_degree, pattern_ris_100(:,9),'-r','LineWidth',2);
grid on;

ylim([0,35]);
xlim([-90,90]);

legend('Radar-only','without RIS','20 elements','60 elements','100 elements','FontSize',12,'interpreter','latex');
xlabel('Azimuth Angle [degree]','FontSize',12,'interpreter','latex');
ylabel('Beampattern Value [dBm]','FontSize',12,'interpreter','latex');












% clc
% clear all
% close all
% addpath('./function/');
% 
% %% trade-off
% load('./ris_aided_data/rician_1000/ris_aided_single_100.mat');
% pattern_ris_100 = pattern_ris;
% wsr_ris = [0; wsr_ris;wsr_ris(end)];
% prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
% plot(wsr_ris, prob_power_ris,'-+r','LineWidth',3);
% hold on;
% 
% 
% load('./ris_aided_data/rician_1000/ris_aided_single_60.mat');
% pattern_ris_60 = pattern_ris;
% wsr_ris = [0; wsr_ris;wsr_ris(end)];
% prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
% plot(wsr_ris, prob_power_ris,'-+c','LineWidth',3);
% 
% load('./ris_aided_data/rician_1000/ris_aided_single.mat');
% pattern_ris_20 = pattern_ris;
% wsr_ris = [0; wsr_ris; wsr_ris(end)];
% prob_power_ris = [prob_power_ris(1); prob_power_ris;0];
% plot(wsr_ris, prob_power_ris,'-+g','LineWidth',3);
% 
% load('baseline.mat');
% wsr_baseline = [0; wsr_baseline;wsr_baseline(end)];
% prob_power_baseline = [prob_power_baseline(1); prob_power_baseline;0];
% plot(wsr_baseline, prob_power_baseline,'-+b','LineWidth',3);
% 
% grid on;
% legend('100 elements','60 elements', '20 elements','without RIS','FontSize',12,'interpreter','latex');
% xlabel('WSR [bps/Hz]','FontSize',12,'interpreter','latex');
% ylabel('Probing Power [dBm]','FontSize',12,'interpreter','latex');
% ylim([19 34]);
% xlim([0 10]);
% 
% %% radar only beampattern
% para = para_init();
% theta_degree = -90:90;
% theta = theta_degree*pi/180;
% pattern_radar_only = zeros(length(theta),1);
% a = ULA_func(0,para.M);
% R = para.Pt/para.M*(a*a');
% 
% for i = 1:length(theta)
%     t = theta(i);
%     a = ULA_func(t,para.M);
%     pattern_radar_only(i) = real(a'*R*a);
% end
% pattern_radar_only = 10*log10(pattern_radar_only);
% 
% 
% %% beampattern WSR = 7.9
% figure;
% plot(theta_degree, pattern_ris_100(:,7),'-r','LineWidth',2);
% hold on;
% plot(theta_degree, pattern_ris_60(:,10),'-c','LineWidth',2);
% plot(theta_degree, pattern_ris_20(:,12),'-g','LineWidth',2);
% plot(theta_degree, pattern_baseline(:,14),'-b','LineWidth',2);
% grid on;
% 
% plot(theta_degree, pattern_radar_only,'--','LineWidth',2,'Color',[0.5 0.5 0.5]);
% ylim([0,35]);
% xlim([-90,90]);
% 
% legend('100 elements','60 elements', '20 elements','without RIS','Radar-only','FontSize',12,'interpreter','latex');
% xlabel('Azimuth Angle [degree]','FontSize',12,'interpreter','latex');
% ylabel('Beampattern Value [dBm]','FontSize',12,'interpreter','latex');

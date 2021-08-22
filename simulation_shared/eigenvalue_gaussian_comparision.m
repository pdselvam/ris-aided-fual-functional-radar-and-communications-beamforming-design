clc
clear all
close all

load('converge_eigen.mat');
plot(1:30, wsr_all,'LineWidth',1);
hold on;

load('converge_gaussian.mat');
plot(1:30, wsr_all,'LineWidth',1);

grid on;
legend('Eigenvalue decomposition','Gaussian randomization','FontSize',12,'interpreter','latex');
ylabel('WSR [bps/Hz]','FontSize',12,'interpreter','latex');
xlabel('Steps','FontSize',12,'interpreter','latex');
xlim([1,30])

ax = axes('Position', [0.6, 0.3, 0.2, 0.1]);
plot(ax, 20:30, wsr_all(20:30),'LineWidth',1, 'Color',[0.8500 0.3250 0.0980]);
grid on;

x = [0.8, 0.75];
y = [0.17, 0.3];
annotation('textarrow', x, y);
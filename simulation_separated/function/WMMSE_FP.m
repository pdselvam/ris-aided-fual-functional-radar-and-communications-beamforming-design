function [P, Rq, Theta, wsr, wsr_all] = WMMSE_FP(para, h, Hc, Hr, dc, dr, ar, Z, rho)
%The WMMSE-FP alternating optimization algorithm
%  [P, Rq, Theta, wsr] = WMMSE_FP(para, h, Hc, Hr, dc, dr, ar, Z, rho)
%Inputs:
%   para: structure of the initial parameters
%   h: RIS-user channels
%   Hc: BS-RIS channel for communication part
%   Hr: BS-RIS channel for radar part
%   dc: BS-user channels for communication part
%   dr: BS-user channels for communication part
%   ar: steering vector of radar antennas
%   Z: a positive semidefinite matrix
%   rho: the regularization parameter
%Outputs:
%   P: optimal P
%   Rq: optimal Rq
%   Theta: optimal Theta
%   wsr: optimal WSR
%   wsr: all wsr in each iteration
%Date: 30/05/2021
%Author: Zhaolin Wang


% Input
P = exp(1i.*rand(para.Mc, para.Mc).*2.*pi);
P = sqrt(para.Pc/para.Mc) * (P./vecnorm(P));

Rq = (para.Pr/para.Mr) * eye(para.Mr);

theta=exp(1i.*rand(para.N,1).*2.*pi);
Theta = diag(theta);

if rho < 1
    epsilon = 0.0001;
else
    epsilon = 0.0005;
end

% Step 1: calculate WSR from P, Rq and Theta
c = Hc' * Theta * h + dc;
r = Hr' * Theta * h + dr;
wsr = WSR(para, c, r, P, Rq);
wsr_all = wsr;
diff = 10;
step = 0;
while diff > epsilon
    wsr_pre =wsr;
    step = step + 1;
    % Step 3: update gk
    g_MMSE = MMSE_receiver(para, c, r, P, Rq);
    
    % Step 4: update wk
    e_MMSE = MSE(para, c, r, P, Rq, g_MMSE);
    w = para.weight ./ e_MMSE; 
    
    % Step 5: update P and Rq with updated wk
    [P,Rq] = update_P_Rq(para, rho, c, r, P, Rq, ar, Z, w, g_MMSE);
    
    % Step 6: Update alpha_k
    alpha = SINR(para, c, r, P, Rq);
    
    % Step 7 & 8: Update y_k and Theta
    Theta = update_y_Theta(para, P, Rq, Theta, h, Hc, Hr, dc, dr, alpha);
    
    c = Hc' * Theta * h + dc;
    r = Hr' * Theta * h + dr;
    wsr = WSR(para, c, r, P, Rq);
    diff = abs(wsr - wsr_pre);
    wsr_all = [wsr_all wsr];
    disp(['step--' num2str(step) ', WSR--' num2str(wsr)]);
end


end


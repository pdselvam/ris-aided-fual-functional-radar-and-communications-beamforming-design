function [P, Theta, wsr] = WMMSE_FP_v2(para, h, H, d, Z, rho)
%The WMMSE-FP alternating optimization algorithm for generalized IRS model
%  [P, Theta, wsr] = WMMSE_FP(para, h, Hc, Hr, dc, dr, ar, Z, rho)
%Inputs:
%   para: structure of the initial parameters
%   h: RIS-user channels
%   H: BS-RIS channel
%   d: BS-user channels
%   Z: a positive semidefinite matrix
%   rho: the regularization parameter
%Outputs:
%   P: optimal P
%   wsr: optimal WSR
%   Theta: optimal Theta
%Date: 14/07/2021
%Author: Zhaolin Wang


% Input
P = exp(1i.*rand(para.M, para.M).*2.*pi);
P = sqrt(para.Pt/para.M) * (P./vecnorm(P));


% Group connected IRS model
Theta = zeros(para.N,para.N);
Ng = para.N / para.G;
for i = 1:para.G
    Theta(Ng*(i-1)+1:Ng*i, Ng*(i-1)+1:Ng*i) = exp(1i.*rand(Ng,Ng).*2.*pi);   
end

if rho < 1
    epsilon = 0.0001;
else
    epsilon = 0.001;
end


% Step 1: calculate WSR from P, Rq and Theta
c = H' * Theta * h + d;
wsr = WSR(para, c, P);
diff = 10;
step = 0;
while diff > epsilon 
    wsr_pre =wsr;
    step = step + 1;
    % Step 3: update gk
    g_MMSE = MMSE_receiver(para, c, P);
    
    % Step 4: update wk
    e_MMSE = MSE(para, c, P, g_MMSE);
    w = para.weight ./ e_MMSE; 
    
    % Step 5: update P and Rq with updated wk
    [P] = update_P(para, rho, c, P, Z, w, g_MMSE);
    
    % Step 6: Update alpha_k
    alpha = SINR(para, c, P);
    
    % Step 7 & 8: Update y_k and Theta
    [Theta] = update_y_Theta_v2(para, P, Theta, h, H, d, alpha);
    
    c = H' * Theta * h + d;
    wsr = WSR(para, c, P);
    diff = abs(wsr - wsr_pre);
end
end


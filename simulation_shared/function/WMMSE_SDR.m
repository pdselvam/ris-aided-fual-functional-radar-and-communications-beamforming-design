function [P, wsr] = WMMSE_SDR(para, d, Z, rho)
%The WMMSE-SDP alternating optimization algorithm
%  [P, Rq, Theta, wsr] = WMMSE_FP(para, h, Hc, Hr, dc, dr, ar, Z, rho)
%Inputs:
%   para: structure of the initial parameters
%   dc: BS-user channels for communication part
%   dr: BS-user channels for communication part
%   ar: steering vector of radar antennas
%   Z: a positive semidefinite matrix
%   rho: the regularization parameter
%Outputs:
%   P: optimal P
%   Rq: optimal Rq
%   wsr: optimal WSR
%Date: 01/06/2021
%Author: Zhaolin Wang

% Input
P = exp(1i.*rand(para.M, para.M).*2.*pi);
P = sqrt(para.Pt/para.M) * (P./vecnorm(P));

if rho < 1
    epsilon = 0.0001;
else
    epsilon = 0.0005;
end


% Step 1: calculate WSR from P, Rq and Theta
wsr = WSR(para, d, P);
diff = 10;
step = 0;
while diff > epsilon 
    wsr_pre =wsr;
    step = step + 1;
    % Step 3: update gk
    g_MMSE = MMSE_receiver(para, d, P);
    
    % Step 4: update wk
    e_MMSE = MSE(para, d, P, g_MMSE);
    w = para.weight ./ e_MMSE; 
    
    % Step 5: update P and Rq with updated wk
    
    [P] = update_P(para, rho, d, P, Z, w, g_MMSE);

    wsr = WSR(para, d, P);
    fprintf('WSR: %f\n', wsr);
    diff = abs(wsr - wsr_pre); 
%     disp(['step--' num2str(step) ', WSR--' num2str(wsr(end))]);
end


end


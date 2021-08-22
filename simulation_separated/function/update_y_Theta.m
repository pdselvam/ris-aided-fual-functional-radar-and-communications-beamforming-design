function [Theta] = update_y_Theta(para, P, Rq, Theta, h, Hc, Hr, dc, dr, alpha)
%Update y and Theta using cvx toolbox
%  [Theta] = update_y_Theta(para, P, Rq, h, Hc, Hr, dc, dr, alpha, Theta)
%Inputs:
%   para: structure of the initial parameters
%   P: beamformer at BS
%   Rq: covariance matrix of radar signal
%   Theta: passive beamformer at RIS
%   h: RIS-user channels
%   Hc: BS-RIS channel for communication part
%   Hr: BS-RIS channel for radar part
%   dc: BS-user channels for communication part
%   dr: BS-user channels for radar part
%   alpha: SINR
%Outputs:
%   Theta: updated Theta
%Date: 31/05/2021
%Author: Zhaolin Wang

% Calculate y, U, and v
U = 0;
v = 0;
theta = diag(Theta);
for k = 1:para.K
    hk = h(:,k);
    dck = dc(:,k);
    drk = dr(:,k);
    alpha_k = alpha(k);
    uk = para.weight(k);

    ak = diag(hk')*Hc*P;
    bk = dck'*P;
    Bk = diag(hk')*Hr*Rq*Hr'*diag(hk);
    fk = diag(hk')*Hr*Rq*drk;
    nk = real(drk'*Rq*drk + para.n);

    % Optimal yk
    yk = ( sqrt(uk*(1+alpha_k)) * (theta'*ak(:,k) + bk(k)) ) ...
        / (sum_term(para.K, theta, ak, bk) + theta'*Bk*theta + 2*real(theta'*fk) + nk);
    
    U = U + abs(yk)^2 * (Bk + ak*ak');
    v = v + conj(yk)*sqrt(uk*(1+alpha_k))*ak(:,k) - abs(yk)^2*(fk + sum(conj(bk).*ak,2));
end

cvx_begin quiet
    variable theta(para.N,1) complex

    % constraint
    abs(theta) <= 1;
    
    % objective function
    [P,C] = eig(U) ;
    U_de = P * (C.^(1/2)) * P';
    h = -theta' * (U_de*U_de') * theta + 2*real(theta'*v);
    
    maximize(h);
cvx_end

Theta = diag(theta);

end

function [sum] = sum_term(K, theta, ak, bk)
    sum = 0;
    for j = 1:K  
        sum = sum + pow_abs(theta'*ak(:,j) + bk(j), 2);   
    end
end


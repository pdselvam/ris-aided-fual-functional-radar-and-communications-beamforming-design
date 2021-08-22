function [Theta] = update_y_Theta(para, P, Theta, h, H, d, alpha)
%Update y and Theta using cvx toolbox
%  [Theta] = update_y_Theta(para, P, Rq, h, Hc, Hr, dc, dr, alpha, Theta)
%Inputs:
%   para: structure of the initial parameters
%   P: beamformer at BS
%   Theta: passive beamformer at RIS
%   h: RIS-user channels
%   H: BS-RIS channel
%   d: BS-user channels
%   alpha: SINR
%Outputs:
%   Theta: updated Theta
%Date: 14/07/2021
%Author: Zhaolin Wang

% Calculate y, U, and v
U = 0;
v = 0;
theta = diag(Theta);
for k = 1:para.K
    hk = h(:,k);
    dk = d(:,k);
    alpha_k = alpha(k);
    uk = para.weight(k);

    ak = diag(hk')*H*P;
    bk = dk'*P;
    nk = para.n;

    % Optimal yk
    yk = ( sqrt(uk*(1+alpha_k)) * (theta'*ak(:,k) + bk(k)) ) ...
        / (sum_term(para.K, theta, ak, bk) + nk);
    
    U = U + abs(yk)^2 * (ak*ak');
    v = v + conj(yk)*sqrt(uk*(1+alpha_k))*ak(:,k) - abs(yk)^2*(sum(conj(bk).*ak,2));
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


function [Theta] = update_y_Theta_v2(para, P, Rq, Theta, h, Hc, Hr, dc, dr, alpha)
%Update y and Theta using Quasi-Newton Method
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
%Date: 15/06/2021
%Author: Zhaolin Wang

% Calculate y, U, and v
U = 0;
v = 0;
theta = Theta.';
theta = theta(:);
for k = 1:para.K
    hk = h(:,k);
    dck = dc(:,k);
    drk = dr(:,k);
    alpha_k = alpha(k);
    uk = para.weight(k);
    
    Ak = zeros(para.N*para.N, para.N);
    hk_ext = [conj(hk); zeros((para.N-1)*para.N,1)];
    for i = 0:para.N-1
        Ak(:,i+1) = circshift(hk_ext, i*para.N);
    end  
%     a = [1:para.N, zeros(1,para.N^2)];
%     a_ = repmat(a, 1, para.N);
%     a_ = a_(1:para.N^3);
%     i = find(reshape(a_', para.N^2, para.N));
%     Ak = zeros(para.N*para.N, para.N);
%     Ak(i) = repmat(hk, 1,20);

    
    ak = Ak*Hc*P;
    bk = dck'*P;
    Bk = Ak*Hr*Rq*Hr'*Ak';
    fk = Ak*Hr*Rq*drk;
    nk = real(drk'*Rq*drk + para.n);

    % Optimal yk
    yk = ( sqrt(uk*(1+alpha_k)) * (theta'*ak(:,k) + bk(k)) ) ...
        / (sum_term(para.K, theta, ak, bk) + theta'*Bk*theta + 2*real(theta'*fk) + nk);
    
    U = U + abs(yk)^2 * (Bk + ak*ak');
    v = v + conj(yk)*sqrt(uk*(1+alpha_k))*ak(:,k) - abs(yk)^2*(fk + sum(conj(bk).*ak,2));
end

Ng = para.N / para.G;
x_init = randn(Ng*(Ng+1)/2 * para.G,1);

% f = obj_func(x_init, para.N, para.G, Ng, U, v);

% options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','MaxIter', 100,'OptimalityTolerance', 1e-5);
options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton','MaxIter', 500,'OptimalityTolerance', 1e-5,'MaxFunctionEvaluations',2.1e+7);
[x] = fminunc(@(x)(obj_func(x, para.N, para.G, Ng, U, v)), x_init, options);

Theta = x2Theta(x,para.N, para.G, Ng);

end

function [f] = obj_func(x, N, G, Ng, U, v)
    
    Theta = x2Theta(x, N, G, Ng);
    theta = Theta.';
    theta = theta(:);
    f = theta' * U * theta - 2*real(theta'*v);
%     f = - 2*real(theta'*v);
    f = real(f);
end

function [Theta] = x2Theta(x, N, G, Ng)
    Z0 = 50;
    Theta = zeros(N,N);
    index=find(tril(ones(Ng)));
    select_matrix = tril(ones(Ng),-1);
    num = Ng*(Ng+1)/2;
    for i = 1:G
        xg = x(num*(i-1)+1:num*i);
        Xg = zeros(Ng);
        Xg(index) = xg; 
        Xg = Xg + (Xg .* select_matrix)';
        Theta_g = inv(1i*Xg + Z0*eye(Ng)) * (1i*Xg - Z0*eye(Ng));
        Theta(Ng*(i-1)+1:Ng*i, Ng*(i-1)+1:Ng*i) = Theta_g;
    end 

end

function [sum] = sum_term(K, theta, ak, bk)
    sum = 0;
    for j = 1:K  
        sum = sum + pow_abs(theta'*ak(:,j) + bk(j), 2);   
    end
end


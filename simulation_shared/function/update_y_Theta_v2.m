function [Theta] = update_y_Theta_v2(para, P, Theta, h, H, d, alpha)
%Update y and Theta using Quasi-Newton Method
%  [Theta] = update_y_Theta(para, P, Rq, h, Hc, Hr, dc, dr, alpha, Theta)
%Inputs:
%   para: structure of the initial parameters
%   P: beamformer at BS
%   Theta: passive beamformer at RIS
%   h: RIS-user channels
%   Hc: BS-RIS channel
%   d: BS-user channels
%   alpha: SINR
%Outputs:
%   Theta: updated Theta
%Date: 14/07/2021
%Author: Zhaolin Wang

% Calculate y, U, and v
U = 0;
v = 0;
theta = Theta.';
theta = theta(:);
for k = 1:para.K
    hk = h(:,k);
    dk = d(:,k);
    alpha_k = alpha(k);
    uk = para.weight(k);
    
    Ak = zeros(para.N*para.N, para.N);
    hk_ext = [conj(hk); zeros((para.N-1)*para.N,1)];
    for i = 0:para.N-1
        Ak(:,i+1) = circshift(hk_ext, i*para.N);
    end 
    
    ak = Ak*H*P;
    bk = dk'*P;
    nk = para.n;

    % Optimal yk
    yk = ( sqrt(uk*(1+alpha_k)) * (theta'*ak(:,k) + bk(k)) ) ...
        / (sum_term(para.K, theta, ak, bk) + nk);
    
    U = U + abs(yk)^2 * (ak*ak');
    v = v + conj(yk)*sqrt(uk*(1+alpha_k))*ak(:,k) - abs(yk)^2*(sum(conj(bk).*ak,2));
end

% Quasi-Newton method
Ng = para.N / para.G;
x_init = randn(Ng*(Ng+1)/2 * para.G,1);
options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton','MaxIter', 500,'OptimalityTolerance', 1e-5,'MaxFunctionEvaluations',2.1e+7);
[x] = fminunc(@(x)(obj_func(x, para.N, para.G, Ng, U, v)), x_init, options);

Theta = x2Theta(x,para.N, para.G, Ng);

end

function [f] = obj_func(x, N, G, Ng, U, v)
    
    Theta = x2Theta(x, N, G, Ng);
    theta = Theta.';
    theta = theta(:);
    f = theta' * U * theta - 2*real(theta'*v);
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


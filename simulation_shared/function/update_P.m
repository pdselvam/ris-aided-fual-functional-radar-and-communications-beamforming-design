function [P,f] = update_P(para, rho, c, P, Z, w, g_MMSE)
%Update P and Rq using cvx toolbox
%  [P,f] = update_P_Rq(para, rho, c, P, Z, w, g_MMSE)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   P: beamformer at BS
%   a: steering vector
%   Z: a positive semidefinite matrix
%   w: weight of MSE
%   g_MMSE: MMSE receiver
%Outputs:
%   P: updated P
%   f: optimal value of objective function
%Date: 14/07/2021
%Author: Zhaolin Wang

cvx_begin quiet
    
    variable T(para.M+1,para.M+1,para.K) complex
    
    % constraint
    for k = 1:para.K
        Tk = T(:,:,k);
        Tk == hermitian_semidefinite(para.M+1);
        Tk(para.M+1, para.M+1) == 1;
    end
    diag(sum(T,3)) == [para.Pt/para.M* ones(para.M,1);para.K];
    
    % objective function
    f = objective_func(para, T, c, Z, w, g_MMSE, rho);
    minimize(real(f));
cvx_end

% % Gaussian randomization
% L = 100; 
% wsr_pre = 0;
% P = zeros(para.M, para.K);
% for l = 1:L
%     P_candidate = randomize_func(para, T);
%     wsr = WSR(para, c, P_candidate);  
%     if wsr > wsr_pre
%         P = P_candidate;
%         wsr_pre = wsr;
%     end
% end

% Eigenvecot approximation
P = eigen_approx(para, T);


end

function [P] = eigen_approx(para, T)

P_ext = zeros(para.M+1,para.K);
for k = 1:para.K
    Tk = T(:,:,k);
    [V,D] = eig(Tk);
    pk = sqrt(D(end,end)) * V(:,end);
    pk = pk / abs(pk(end)); % normalize the last entry to 1
    P_ext(:,k) = pk;
end

power = diag(P_ext(1:end-1,:)*P_ext(1:end-1,:)');
scale = power / (para.Pt/para.M);
P_ext(1:end-1,:) = P_ext(1:end-1,:) ./ sqrt(scale);

P = P_ext(1:end-1,:) ./ P_ext(end,:);
end

function [P] = randomize_func(para, T)

P_ext = zeros(para.M+1,para.K);
for k = 1:para.K
    Tk = T(:,:,k);
    pk = mvnrnd(zeros(para.M+1,1), Tk).';
    pk = pk / abs(pk(end)); % normalize the last entry to 1
    P_ext(:,k) = pk;
end

power = diag(P_ext(1:end-1,:)*P_ext(1:end-1,:)');
scale = power / (para.Pt/para.M);
P_ext(1:end-1,:) = P_ext(1:end-1,:) ./ sqrt(scale);

P = P_ext(1:end-1,:) ./ P_ext(end,:);

end

function [f] = objective_func(para, T, c, Z, w, g_MMSE, rho)

f = 0;
Z = [Z zeros(para.M,1); zeros(para.M,1)' 0];
for k = 1:para.K
    Tk = T(:,:,k);
    wk = w(k);
    ck = c(:,k);
    gk = g_MMSE(k);
    Ck1 = [abs(gk)^2*(ck*ck') zeros(para.M,1); zeros(para.M,1)' 0];
    Ck2 = [abs(gk)^2*(ck*ck') -conj(gk)*ck; -gk*ck' 0];
    
    term_1 = 0;
    for j = 1:para.K
        if j ~= k
           Tj = T(:,:,j);
           term_1 = term_1 + trace(Ck1*Tj);            
        end
    end
    term_2 = trace(Ck2*Tk);
    term_3 = trace(Z*Tk);
    
%     f = f + sqrt(rho/(1+rho))*wk*(term_1 + term_2) + sqrt(1/(1+rho))*term_3;
    f = f + rho*wk*(term_1 + term_2) + term_3;   
end
    
end



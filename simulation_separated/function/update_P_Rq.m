function [P,Rq,f] = update_P_Rq(para, rho, c, r, P, Rq, ar, Z, w, g_MMSE)
%Update P and Rq using cvx toolbox
%  [P,Rq, f] = update_P_Rq(para, rho, c, r, P, Rq, ar, Z, w, g_MMSE)
%Inputs:
%   para: structure of the initial parameters
%   c: equivalent communication channel
%   r: equivalent radar channel
%   P: beamformer at BS
%   Rq: covariance matrix of radar signal
%   ar: steering vector of radar antennas
%   Z: a positive semidefinite matrix
%   w: weight of MSE
%   g_MMSE: MMSE receiver
%Outputs:
%   P: updated P
%   Rq: updated Q
%   f: optimal value of objective function
%Date: 30/05/2021
%Author: Zhaolin Wang

cvx_begin quiet
    variable P(para.Mc,para.K) complex
    variable Rq(para.Mr,para.Mr) hermitian semidefinite complex

    % constraint
    diag(Rq) == para.Pr/para.Mr;
    quad_over_lin(vec(P),1) <= para.Pc;
    

    % objective function
    f = objective_func(para, w, c, r, P, Rq, g_MMSE, ar, Z, rho);

    minimize(f);
cvx_end
f = objective_func(para, w, c, r, P, Rq, g_MMSE, ar, Z, rho);
end

function [f] = objective_func(para, w, c, r, P, Rq, g_MMSE, ar, Z, rho)
    term_1 = w'*MSE(para, c, r, P, Rq, g_MMSE);
    term_2 = - ar'*Rq*ar;
    term_3 = 0;
    for k = 1:para.K
        pk = P(:,k);
        term_3 = term_3 + pk'* Z * pk;
    end

    f = rho * term_1 +  (term_2 + term_3);
end


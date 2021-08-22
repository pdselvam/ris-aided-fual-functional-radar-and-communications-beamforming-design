function [h, H, d, LOS, NLOS] = generate_channel(para, angle, path_loss, a)
%Generate the BS-user, RIS-user and BS-RIS channels 
%  [h, Hc, Hr, dc, dr, LOS, NLOS] = generate_channel(para, angle, path_loss, ar, ac)
%Inputs:
%   para: structure of the initial parameters
%   angle: struture of the angles
%   path_loss: structure of the path loss
%   a: steering vector
%Outputs:
%   h: RIS-user channels
%   H: BS-RIS channel
%   d: BS-user channels
%   LOS: structure of the LOS components of the Rician channel
%   NLOS: structure of the NLOS components of the Rician channel
%Date: 14/07/2021
%Author: Zhaolin Wang

%% Rayleigh channel (BS-user channel)
d = path_loss.BU .* 1/sqrt(2) .* ( randn(para.M,para.K) + 1i*randn(para.M,para.K) );

%% Rican channel (BS-RIS and RIS-user channel)
epsilon = para.rician; % Rician factor

% NLOS components
H_NLOS = 1/sqrt(2) .* ( randn(para.N, para.M) + 1i*randn(para.N, para.M) );
h_NLOS = 1/sqrt(2) .* ( randn(para.N,para.K) + 1i*randn(para.N,para.K) );

a_ris_bs = ULA_func(angle.BS, para.N);
H_LOS = a_ris_bs*a';
H = (sqrt(epsilon/(epsilon+1)) * H_LOS + sqrt(1/(epsilon+1)) * H_NLOS);

h_LOS = [];
for k = 1:para.K
    h_LOS = [h_LOS ULA_func(angle.user(k), para.N)];
       
end
h = path_loss.BRU .* (sqrt(epsilon/(epsilon+1)) * h_LOS + sqrt(1/(epsilon+1)) * h_NLOS);

LOS = struct('H', H_LOS,'h', h_LOS);
NLOS = struct('H', H_NLOS,'h', h_NLOS);

end


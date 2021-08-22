function [h, Hc, Hr, dc, dr, NLOS] = update_channel(para, LOS, NLOS, dc, dr, path_loss)
%Generate the BS-user, RIS-user and BS-RIS channels 
%  [h, Hc, Hr, dc, dr, NLOS] = update_channel(para, LOS, NLOS, dc, dr, path_loss)
%Inputs:
%   para: structure of the initial parameters
%   LOS: structure of the LOS components of the Rician channel
%   NLOS: structure of the NLOS components of the Rician channel
%   dc: BS-user channels for communication part
%   dr: BS-user channels for communication part
%   path_loss: structure of the path loss
%Outputs:
%   h: RIS-user channels
%   Hc: BS-RIS channel for communication part
%   Hr: BS-RIS channel for radar part
%   dc: BS-user channels for communication part
%   dr: BS-user channels for communication part
%   NLOS: structure of the NLOS components of the Rician channel
%Date: 30/05/2021
%Author: Zhaolin Wang

rho = 0.7;
epsilon = para.rician;

%% Rayleigh channel (BS-user channel)
dc_new = path_loss.BU .* 1/sqrt(2) .* ( randn(para.Mc,para.K) + 1i*randn(para.Mc,para.K) );
dr_new = path_loss.BU .* 1/sqrt(2) .* ( randn(para.Mr,para.K) + 1i*randn(para.Mr,para.K) );
dc = rho * dc + sqrt(1-rho^2) * dc_new;
dr = rho * dr + sqrt(1-rho^2) * dr_new;

%% Rician channel

% NLOS components
Hc_NLOS_new = 1/sqrt(2) .* ( randn(para.N, para.Mc) + 1i*randn(para.N, para.Mc) );
Hr_NLOS_new = 1/sqrt(2) .* ( randn(para.N, para.Mr) + 1i*randn(para.N, para.Mr) );
h_NLOS_new = 1/sqrt(2) .* ( randn(para.N,para.K) + 1i*randn(para.N,para.K) );

Hc_NLOS = rho * NLOS.Hc + sqrt(1-rho^2) * Hc_NLOS_new;
Hr_NLOS = rho * NLOS.Hr + sqrt(1-rho^2) * Hr_NLOS_new;
h_NLOS = rho * NLOS.h + sqrt(1-rho^2) * h_NLOS_new;

Hc = (sqrt(epsilon/(epsilon+1)) * LOS.Hc + sqrt(1/(epsilon+1)) * Hc_NLOS);
Hr = (sqrt(epsilon/(epsilon+1)) * LOS.Hr + sqrt(1/(epsilon+1)) * Hr_NLOS);
h = path_loss.BRU .* (sqrt(epsilon/(epsilon+1)) * LOS.h + sqrt(1/(epsilon+1)) * h_NLOS);

NLOS = struct('Hc', Hc_NLOS,'Hr', Hr_NLOS,'h', h_NLOS);

end


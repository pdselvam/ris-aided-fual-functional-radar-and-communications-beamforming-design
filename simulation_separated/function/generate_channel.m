function [h, Hc, Hr, dc, dr, LOS, NLOS] = generate_channel(para, angle, path_loss, ar, ac)
%Generate the BS-user, RIS-user and BS-RIS channels 
%  [h, Hc, Hr, dc, dr, LOS, NLOS] = generate_channel(para, angle, path_loss, ar, ac)
%Inputs:
%   para: structure of the initial parameters
%   angle: struture of the angles
%   path_loss: structure of the path loss
%   ar: steering vector of the radar antennas
%   ac: steering vector of the communication antennas
%Outputs:
%   h: RIS-user channels
%   Hc: BS-RIS channel for communication part
%   Hr: BS-RIS channel for radar part
%   dc: BS-user channels for communication part
%   dr: BS-user channels for radar part
%   LOS: structure of the LOS components of the Rician channel
%   NLOS: structure of the NLOS components of the Rician channel
%Date: 30/05/2021
%Author: Zhaolin Wang

%% Rayleigh channel (BS-user channel)
dc = path_loss.BU .* 1/sqrt(2) .* ( randn(para.Mc,para.K) + 1i*randn(para.Mc,para.K) );
dr = path_loss.BU .* 1/sqrt(2) .* ( randn(para.Mr,para.K) + 1i*randn(para.Mr,para.K) );


%% Rican channel (BS-RIS and RIS-user channel)
epsilon = para.rician; % Rician factor

% NLOS components
Hc_NLOS = 1/sqrt(2) .* ( randn(para.N, para.Mc) + 1i*randn(para.N, para.Mc) );
Hr_NLOS = 1/sqrt(2) .* ( randn(para.N, para.Mr) + 1i*randn(para.N, para.Mr) );
h_NLOS = 1/sqrt(2) .* ( randn(para.N,para.K) + 1i*randn(para.N,para.K) );

a_ris_bs = ULA_func(angle.BS, para.N);
Hc_LOS = a_ris_bs*ac';
Hr_LOS = a_ris_bs*ar';
Hc = (sqrt(epsilon/(epsilon+1)) * Hc_LOS + sqrt(1/(epsilon+1)) * Hc_NLOS);
Hr = (sqrt(epsilon/(epsilon+1)) * Hr_LOS + sqrt(1/(epsilon+1)) * Hr_NLOS);    

h_LOS = [];
for k = 1:para.K
    h_LOS = [h_LOS ULA_func(angle.user(k), para.N)];
       
end
h = path_loss.BRU .* (sqrt(epsilon/(epsilon+1)) * h_LOS + sqrt(1/(epsilon+1)) * h_NLOS);

LOS = struct('Hc', Hc_LOS,'Hr', Hr_LOS,'h', h_LOS);
NLOS = struct('Hc', Hc_NLOS,'Hr', Hr_NLOS,'h', h_NLOS);


end


function [chamber_L, contract_L, nozzle_L, cham_chan_num] = getChamberSize(A_t, contrac, expand, L_star, chan_ID, chan_t);

%% Calculating Chamber Length
A_c = contrac * A_t;
V_c = L_star * A_t;

chamber_L = V_c/A_c;

%% Calculating Contacting Length
%Assumes simple 45 degree contracting angle

r_c = sqrt(A_c/pi);
r_t = sqrt(A_t/pi);

contract_L = r_c - r_t;

%%  Calculating Nozzle Length
%Assumes simple 15 degree conical nozzle

A_e = A_t * expand;

r_e = sqrt(A_e/pi);

nozzle_L = (r_e - r_t)/tand(15);

%% Chamber Channel Number
chan_OD = chan_ID + 2*chan_t;

new_diam = 2*r_c + chan_OD/2;
new_circum = new_diam * pi;
cham_chan_num = new_circum/chan_OD;

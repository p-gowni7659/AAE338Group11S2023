function [chamber_L, contract_L, nozzle_L] = getChamberSize(A_t, contrac, expand, L_star);

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
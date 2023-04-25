function [hgas,area,Tgas,tubelen] = nozzleprops(Dt, Area_arr, h_g_x, wallt, tubenum,T_gas_arr,converge_num,contract_L,diverge_num,nozzle_L)

%Find Area at a distance x from the start of the converging section
x = (tubenum*Dt) - (Dt/2);
% lenarrconv = linspace(0,contract_L, converge_num);
% lenarrdiv = linspace(0,nozzle_L, diverge_num);
% lenarr = [lenarrconv lenarrdiv];
disp(floor(converge_num * x / contract_L))
if x <= contract_L
    area = Area_arr(floor(converge_num * x / contract_L));
    Tgas = T_gas_arr(floor(converge_num * x / contract_L));
    hgas = h_g_x(floor(converge_num * x / contract_L));
elseif x <= contract_L + nozzle_L
    area = Area_arr(floor(diverge_num * x / nozzle_L)+converge_num);
    Tgas = T_gas_arr(floor(diverge_num * x / nozzle_L)+converge_num);
    hgas = h_g_x(floor(diverge_num * x / nozzle_L)+converge_num);
else
    disp("Spence Sucks")
end
nozzlerad = sqrt(area/pi);
tubelen = 2*pi*(nozzlerad + wallt + Dt/2);
disp(hgas)
disp(area)
disp(Tgas)
disp(tubelen)
end


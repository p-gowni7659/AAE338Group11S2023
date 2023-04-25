function [hgas,area,Tgas,tubelen] = nozzleprops(Dt, Area_arr, h_g_x, wallt, tubenum,T_gas_arr,converge_num,contract_L,diverge_num,nozzle_L)

%Find Area at a distance x from the start of the converging section
x = (tubenum*Dt) - (Dt/2);

disp(floor(converge_num * x / contract_L))
run = 1;
if  floor(diverge_num * x / nozzle_L)+converge_num > length(Area_arr)
    area = Area_arr(end);
    Tgas = T_gas_arr(end);
    hgas = h_g_x(end);
    run = 0;
end

if x <= contract_L && run == 1
    area = Area_arr(floor(converge_num * x / contract_L));
    Tgas = T_gas_arr(floor(converge_num * x / contract_L));
    hgas = h_g_x(floor(converge_num * x / contract_L));
elseif x <= contract_L + nozzle_L && run == 1
    area = Area_arr(floor(diverge_num * x / nozzle_L)+converge_num);
    Tgas = T_gas_arr(floor(diverge_num * x / nozzle_L)+converge_num);
    hgas = h_g_x(floor(diverge_num * x / nozzle_L)+converge_num);
elseif run == 1
    disp("Spence Sucks")
end
nozzlerad = sqrt(area/pi);
tubelen = 2*pi*(nozzlerad + wallt + Dt/2);
disp(hgas)
disp(area)
disp(Tgas)
disp(tubelen)
end


function [] = nozzleprops(Dt, Area_arr, gma, M_x, tubenum,T_gas)

%Find Area at a distance x from the start of the converging section
x = (tubenum*Dt) - (Dt/2);

if x <= contract_L
    area = Area_arr(converge_num * x / contract_L);
elseif x <= contract_L + nozzle_L
    area = Area_arr(diverge_num * x / nozzle_L);
else
    disp("Spence Sucks")
end
Area_arr = 




end

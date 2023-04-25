function [] = nozzleprops(Dt, Area_arr, gma, Mach_arr, tubenum)

%Find Area at a distance x from the start of the converging section
x = tubenum*Dt;

if x <= contract_L
    area = Area_arr(converge_num * x / contract_L);
elseif x <= contract_L + nozzle_L
    area = Area_arr(diverge_num * x / nozzle_L);
else
    disp("Spence Sucks")
end



Pr = (4 * gma) / ((9 * gma) - 5);
r = Pr ^ (0.33);

T_gas = T_cns .* (1 + (r .* ((gma - 1) ./ 2) .* M_x)) ./ (1 + (((gma - 1) ./ 2) .* M_x));


end

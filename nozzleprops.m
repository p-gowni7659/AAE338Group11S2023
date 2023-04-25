function [] = nozzleprops(Dt, AreaRatio_arr, At, gma, Mach_arr, tubenum)


AreaRatio_arr = 

Pr = (4 * gma) / (9 * gma - 5);
r = Pr ^ (0.33);

T_gas = T_cns .* (1 + (r .* ((gma - 1) ./ 2) .* M_x)) ./ (1 + (((gma - 1) ./ 2) .* M_x));



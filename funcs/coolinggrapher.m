function [] = coolinggrapher(M_hel_arr, qdot_arr, T_hw_arr, T_cw_arr, T_hel_arr, P_hel_arr, Aratio, x, radius)

figure()
plot(M_hel_arr, T_hel_arr)
grid on
xlabel("Mach Number")
ylabel("Helium Temperature")

figure()
plot(M_hel_arr, T_hw_arr)
grid on
xlabel("Mach Number")
ylabel("Chamber Side Wall Temperature")

figure()
plot(M_hel_arr, T_cw_arr)
grid on
xlabel("Mach Number")
ylabel("Far Side Wall Temperature")

figure()
plot(M_hel_arr, qdot_arr)
grid on
xlabel("Mach Number")
ylabel("Specific heat transfer")

figure()
plot(M_hel_arr, P_hel_arr)
grid on
xlabel("Mach Number")
ylabel("Helium Static Pressure")




% Plot the cross-section
figure()
plot(xplot, sqrt(A/pi));
grid on
hold on
plot(x3, x4, 'color', [0, 0.4470, 0.7410])
yline(0)
xlabel('Radial distance (m)');
ylabel('y (m)');
title('Rocket Nozzle and Chamber Cross-Section');
axis equal


figure()
plot(Aratio, T_gas)
grid on;
title('Nozzle Temperature Gas Distribution [Correlation]')
xlabel("Area Ratio [A/At]")
ylabel('Gas Temperature [K]')

figure()
plot(Aratio, h_g_x)
grid on;
title('Nozzle Convective Heat Transfer Coeff Distribution')
xlabel("Area Ratio [A/At]")

figure()
plot(Aratio, Qdot_x)
grid on;
title('Nozzle Qdot Distribution')
xlabel("Area Ratio [A/At]")
ylabel('Qdot')




plot(Aratio, M_x)
grid on;
title('Nozzle Mach Number Distribution')
xlabel("Area Ratio [A/At]")
ylabel('Mach Number')

figure()
plot(Aratio, rho_x)
grid on;
title('Nozzle Density Distribution')
xlabel("Area Ratio [A/At]")
ylabel('Density [kg/m^3]')

figure()
plot(Aratio, T_x)
grid on;
title('Nozzle Temperature Distribution [Isentropic]')
xlabel("Area Ratio [A/At]")
ylabel('Temperature [K]')

figure()
plot(Aratio, V_x)
grid on;
title('Nozzle Velocity Distribution')
xlabel("Area Ratio [A/At]")
ylabel('Velocity [m/s]')
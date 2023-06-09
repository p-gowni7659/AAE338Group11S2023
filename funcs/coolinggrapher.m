function [] = coolinggrapher(M_hel_arr, qdot_arr, T_hw_arr, T_cw_arr, T_hel_arr, P_hel_arr, Aratio, xplot, A, T_gas, h_g_x, Qdot_x, M_x, rho_x, T_x, V_x, x3, x4, T_hw_arr_cha, T_hw, M_hel_arr_cha, P_hel_arr_cha, qdot_arr_cha, rho_arr_cha, T_hw_arr_nozz, rho_arr, h_gas)

figure()
plot(linspace(0,100,length(M_hel_arr)), M_hel_arr)
grid on
xlabel("Distance Normalized to Loop Length (%)")
ylabel("Mach Number")
title('Helium Mach Number Distribution in Chamber Loop')

figure()
plot(linspace(0,100,length(rho_arr)), rho_arr)
grid on
xlabel("Distance Normalized to Loop Length (%)")
ylabel("Helium Density [kg/m^3]")
title('Helium Density Distribution in Chamber Loop')

plot(linspace(0,100,length(T_hel_arr)), T_hel_arr)
grid on
xlabel("Distance Normalized to Loop Length (%)")
ylabel("Helium Temperature [K]")
title('Helium Temperature Distribution in Chamber Loop')

figure()
plot(linspace(0,100,length(T_hel_arr)), T_hw_arr)
grid on
xlabel("Distance Normalized to Loop Length (%)")
ylabel("Chamber Side Wall Temperature [K]")
title('Chamber Wall Temperature Distribution in Chamber Loop')

figure()
plot(linspace(0,100,length(T_hel_arr)), T_cw_arr)
grid on
xlabel("Distance Normalized to Loop Length (%)")
ylabel("Channel Wall Temperature [K]")
title('Channel Wall Temperature Distribution in Chamber Loop')

figure()
plot(linspace(0,100,length(T_hel_arr)), qdot_arr)
grid on
xlabel("Distance Normalized to Loop Length (%)")
ylabel("Heat Flux [W/(m^2)]")
title("Steady State Heat Flux Distribution in Chamber Loop")

figure()
plot(linspace(0,100,length(T_hel_arr)), P_hel_arr./1000)
grid on
xlabel("Distance Normalized to Loop Length (%)")
ylabel("Helium Static Pressure [kPa]")
title('Helium Static Pressure Distribution in Chamber Loop')




% Plot the cross-section
figure()
plot(xplot, sqrt(A/pi));
grid on
hold on
plot(x3, x4, 'color', [0, 0.4470, 0.7410])
yline(0)
xlabel('Linear Distance Referenced from Nozzle Throat [m]');
ylabel('Radial Distance [m]');
title('Rocket Nozzle and Chamber Cross-Section');
axis equal


% figure()
% plot(Aratio, T_gas)
% grid on;
% title('Nozzle Temperature Gas Distribution [Correlation]')
% xlabel("Area Ratio [A/At]")
% ylabel('Gas Temperature [K]')
% 
% figure()
% plot(Aratio, h_g_x)
% grid on;
% title('Nozzle Convective Heat Transfer Coeff Distribution')
% xlabel("Area Ratio [A/At]")
% 
% figure()
% plot(Aratio, Qdot_x)
% grid on;
% title('Nozzle Qdot Distribution')
% xlabel("Area Ratio [A/At]")
% ylabel('Qdot')
% 
% 
% 
% 
% plot(Aratio, M_x)
% grid on;
% title('Nozzle Mach Number Distribution')
% xlabel("Area Ratio [A/At]")
% ylabel('Mach Number')
% 
% figure()
% plot(Aratio, rho_x)
% grid on;
% title('Nozzle Density Distribution')
% xlabel("Area Ratio [A/At]")
% ylabel('Density [kg/m^3]')
% 
% figure()
% plot(Aratio, T_x)
% grid on;
% title('Nozzle Temperature Distribution [Isentropic]')
% xlabel("Area Ratio [A/At]")
% ylabel('Temperature [K]')
% 
% figure()
% plot(Aratio, V_x)
% grid on;
% title('Nozzle Velocity Distribution')
% xlabel("Area Ratio [A/At]")
% ylabel('Velocity [m/s]')

figure()
find_zero = find(T_hw_arr_cha(1,:) == 0);
plot(linspace(0,100,length(T_hw_arr_cha(1,1:(find_zero(1) - 1)))), T_hw_arr_cha(1,1:(find_zero(1) - 1)))
grid on
hold on
i = 1;
while i < size(T_hw_arr_cha,1)
    find_zero = find(T_hw_arr_cha(i+1,:) == 0);
    if isempty(find_zero)
        find_zero = size(T_hw_arr_cha,2);
    end
    plot(linspace(0,100,length(T_hw_arr_cha(i+1,1:(find_zero(1) - 1)))), T_hw_arr_cha(i+1,1:(find_zero(1) - 1)))
    i = i+1;
end
yline(1088,"LineWidth",3)
xlabel("Normalized Channel Length")
ylabel("Chamber Side Temperature [K]")
title("Chamber Side Temperature vs Channel Length in Nozzle", "For Channels in Converging/Diverging Section")


find_zero_new = ones(10000,1);
for i = 1:size(M_hel_arr_cha,1)
    find_zero = find(isnan(M_hel_arr_cha(i,:)));
    if length(find_zero) > length(find_zero_new)
        find_zero_final = find_zero;
        i_final = i;
    end
    find_zero_new = find_zero;
end

figure()
plot(linspace(0,100,length(M_hel_arr_cha(i_final,1:(find_zero_final(1) - 1)))),M_hel_arr_cha(i_final,1:(find_zero_final(1) - 1)))
grid on
xlabel("Normalized Channel Length (%)")
ylabel("Channel Mach Number")
title("Mach Number vs Channel Length for Throat")

figure()
plot(linspace(0,100,length(P_hel_arr_cha(i_final,1:(find_zero_final(1) - 1)))),P_hel_arr_cha(i_final,1:(find_zero_final(1) - 1)) / 1000)
grid on
xlabel("Normalized Channel Length (%)")
ylabel("Channel Pressure (KPa")
title("Pressure vs Channel Length for Throat")

figure()
plot(linspace(0,100,length(qdot_arr_cha(i_final,1:(find_zero_final(1) - 1)))),qdot_arr_cha(i_final,1:(find_zero_final(1) - 1)))
grid on
xlabel("Normalized Channel Length (%)")
ylabel("Channel Heat Flux (W/m^2)")
title("Heat Flux vs Channel Length for Throat")

figure()
plot(linspace(0,100,length(rho_arr_cha(i_final,1:(find_zero_final(1) - 1)))),rho_arr_cha(i_final,1:(find_zero_final(1) - 1)))
grid on
xlabel("Normalized Channel Length (%)")
ylabel("Channel Density (kg/m^3)")
title("Density vs Channel Length for Throat")

figure()
plot(linspace(0,100,length(T_hw_arr_cha(i_final,1:(find_zero_final(1) - 1)))),T_hw_arr_cha(i_final,1:(find_zero_final(1) - 1)))
grid on
xlabel("Normalized Channel Length (%)")
ylabel("Channel Temperature (K)")
title("Temperature vs Channel Length for Throat")



function [] = coolinggrapher(M_hel_arr,qdot_arr,T_hw_arr,T_cw_arr,T_hel_arr, P_hel_arr)

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
plot(M_hel_arr, qdot_arr)
grid on
xlabel("Mach Number")
ylabel("Specific heat transfer")

figure()
plot(M_hel_arr, P_hel_arr)
grid on
xlabel("Mach Number")
ylabel("Helium Static Pressure")
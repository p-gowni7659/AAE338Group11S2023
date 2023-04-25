%AAE338 Project

close all;
clear;
clc;

%% Initialization
pressureExit = 2; %psia
pressureChamber = 300; %psia
OF = 2.35;
mdot = 5; %kg/s propellant mass flow rate (iterate this variable?)
contractionRatio = 3;
Lstar = 1.143; %meters, 45 in
T_hw = 1473.15; %K -Inconcel X750 Wall Temperature
k = 25; %W/m-K -Nozzle wall made from Inconel X750

%% CEA
CEAPath = append(pwd, '/PSP_CEA_function_wrapper');
INPPath = append(pwd, '/INP_OUT');
funcPath = append(pwd, '/funcs');
addpath(CEAPath);
addpath(INPPath);
addpath(funcPath);

nameString = strcat('338_estimates_pip_', num2str(int8(pressureChamber / pressureExit)), '_p_c_', num2str(pressureChamber), '_O_F_', num2str(OF));

inputName = append(nameString, '.inp');
outputName = append(nameString, '.out');

[Isp, CStar, expansionRatio, specificHeatRatio, combustionTemperature, cp, conductivity, enthalpy, rho0] = PSP_1DOF_CEA_function_wrapper(pressureChamber,pressureExit, OF, nameString, 0);
movefile(inputName, 'INP_OUT');
movefile(outputName, 'INP_OUT');
delete(append(pwd, '\PSP_CEA_function_wrapper\', inputName));

Isp = Isp / 9.81; %seconds

gma = specificHeatRatio;
P0 = pressureChamber * 6894.76;
T_cns = combustionTemperature; %K

R = P0 / (rho0 * T_cns);

Aratio_sub = linspace(contractionRatio, 1, 200);
Aratio_sup = linspace(1.001, expansionRatio, 3000);
Aratio = [Aratio_sub Aratio_sup];
M_x = [];
rho_x = [];
T_x = [];
for subratio = Aratio_sub
    [mach_sub, T_sub, P_sub, rho_sub, area_sub] = flowisentropic(gma, subratio, 'sub');
    M_x = [M_x mach_sub];
    rho_x = [rho_x rho_sub];
    T_x = [T_x T_sub];
end
for supratio = Aratio_sup
    [mach_sup, T_sup, P_sup, rho_sup, area_sup] = flowisentropic(gma, supratio, 'sup');
    M_x = [M_x mach_sup];
    rho_x = [rho_x rho_sup];
    T_x = [T_x T_sup];
end
rho_x = rho_x .* rho0;
T_x = T_x .* T_cns;
V_x = M_x .* sqrt(gma .* R .* T_x);

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
%%
Pr = (4 * gma) / (9 * gma - 5);
r = Pr ^ (0.33);

T_gas = T_cns .* (1 + (r .* ((gma - 1) ./ 2) .* M_x)) ./ (1 + (((gma - 1) ./ 2) .* M_x));

h_g_x = (rho_x .* V_x) .^ 0.8;
Qdot_x = h_g_x .* (T_gas - T_hw);
%This Qdot_x must be matched by the regenerative cooling from gas(rayliegh flow)

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

%% Chamber

chan_ID = 0;
chan_t = 0;
At = (mdot * CStar) / P0; %m^2

[chamber_L, contract_L, nozzle_L, cham_chan_num] = getChamberSize(At, contractionRatio, expansionRatio, Lstar, chan_ID, chan_t);

Dt = 2 * sqrt(At/pi);
Dc = sqrt((At * contractionRatio) / pi) * 2;
De = sqrt((At * expansionRatio) / pi) * 2;

A = ((At ./ M_x) .* (((2 + (gma - 1) .* M_x .^ 2) ./ (gma + 1)) .^ ((gma + 1) ./ (2 .* (gma - 1)))));

x1 = linspace(-contract_L, 0, 200);
x2 = linspace(0, nozzle_L, 3000);
x = [x1 x2];
x3 = linspace(-chamber_L, -contract_L, 5);
x4 = sqrt(A(1)/pi) * ones(length(x3));

% Plot the cross-section
figure()
plot(x, sqrt(A/pi));
grid on
hold on
plot(x3, x4, 'color', [0, 0.4470, 0.7410])
yline(0)
xlabel('Radial distance (m)');
ylabel('y (m)');
title('Rocket Nozzle and Chamber Cross-Section');
axis equal

%% Helium Initial Conditions/Loop
%Initial Helium Pressure
Chambertemp = 3000; %K

heltemp_init = 200; %K
helpress_init = 5000000;%Pa
helmach_init = 0.01;
h_gas = 100;
k_wall = 25;
wall_thick = 0.02;
Dh = 0.01;
Tstag_hel_init = 3487;
%% LOOP
step = 0.01;
i = 1;

qdot_arr = [];
T_cw_arr = [];
T_hw_arr = [];
T_hel_arr = [heltemp_init];
P_hel_arr = [helpress_init];
M_hel_arr = [helmach_init];

while i < (chamber_L/step)

    %Setting Helium Step Input Parameters
    if i == 1
        Ti = heltemp_init;
        Pi = helpress_init;
        Mi = helmach_init;
        T0i = Tstag_hel_init;
    else
        Ti = Te;
        Pi = Pe;
        Mi = Me;
        T0i = T0e;
    end
    T_hel_arr = [T_hel_arr Ti];
    P_hel_arr = [P_hel_arr Pi];
    M_hel_arr = [M_hel_arr Mi];
    %Initializing Forced Convection
    Cp = py.CoolProp.CoolProp.PropsSI("C","T",Ti,"P", Pi,"Helium");
    Cv = py.CoolProp.CoolProp.PropsSI("O","T",Ti,"P", Pi,"Helium");
    gma_hel = Cp/Cv;

    R_i = py.CoolProp.CoolProp.PropsSI("gas_constant","T",Ti,"P", Pi,"Helium");
    Vel_i = Mi * sqrt(gma_hel * R_i * Ti);
    [q_dot, T_cw, T_hw] = convergeTemp(Chambertemp, h_gas, k, wall_thick, Ti,Vel_i,Dh,Pi);

    %Checking Forced Convection Convergence
    if q_dot == 1
       
        disp('Forced Convection Failed to Converge')
        break
    end

    %Storing Forced Convection Results
    qdot_arr = [qdot_arr q_dot];
    T_cw_arr = [T_cw_arr T_cw];
    T_hw_arr = [T_hw_arr T_hw];
    
    %Initialiing Ralyeigh Flow
    Qdot = q_dot * step * (Dh+2*wall_thick);
    rho_i = py.CoolProp.CoolProp.PropsSI("D","T",Ti,"P", Pi,"Helium");
    mdot = rho_i * Vel_i * (pi*(Dh/2)^2);
    
    [T0e] = getTempStagNew(Qdot, mdot, Mi, Ti, Pi, gma_hel, Cp);
    [Me,Te, Pe] = RayleighFlow(Pi, Ti, T0e, Mi, gma_hel);
end

%% Bottom of Script
% Resets Matlabs Path preference so it doesn't mess up your matlab
clear CEApath INPPath funcPath

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED
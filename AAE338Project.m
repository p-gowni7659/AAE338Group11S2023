%AAE338 Project

close all;
clear;
clc;


% Using CoolProp
py.CoolProp.CoolProp.PropsSI('P', 'T', 298, 'Q', 0, 'water');


%% CEA
CEAPath = append(pwd, '/PSP_CEA_function_wrapper');
INPPath = append(pwd, '/INP_OUT');
funcPath = append(pwd, '/funcs');
addpath(CEAPath);
addpath(INPPath);
addpath(funcPath);

pressureExit = 5; %psia
pressureChamber = 150; %psia
OF = 2.35;
mdot = 5; %kg/s propellant mass flow rate (iterate this variable?)
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
%As an example, used F1 Engine Epxanison Ratio of 16:1. Guessed Chamber
%Area Ratio of 3:1 since it will change.
Aratio_sub = linspace(3,1,200);
Aratio_sup = linspace(1.01,expansionRatio,3000);
Aratio = [Aratio_sub,Aratio_sup];
M_x = [];
rho_x = [];
T_x = [];
for subratio = Aratio_sub
    [mach_sub,T_sub,P_sub,rho_sub,area_sub] = flowisentropic(gma,subratio,'sub');
    M_x = [M_x mach_sub];
    rho_x = [rho_x rho_sub];
    T_x = [T_x T_sub];
end
for supratio = Aratio_sup
    [mach_sup,T_sup,P_sup,rho_sup,area_sup] = flowisentropic(gma,supratio,'sup');
    M_x = [M_x mach_sup];
    rho_x = [rho_x rho_sup];
    T_x = [T_x T_sup];
end
rho_x = rho_x .* rho0;
T_x = T_x .* T_cns;
V_x = M_x .* sqrt(gma.*R.*T_x);
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
Pr = (4*gma)/(9*gma - 5);
r = Pr^(0.33);

T_gas = T_cns .* (1+(r.*((gma-1)./2).*M_x))./(1+(((gma-1)./2).*M_x));

h_g_x = (rho_x .* V_x).^0.8;
T_hw = 1473.15; %K -Inconcel X750 Wall Temperature
k = 25; %W/m-K -Nozzle wall made from Inconel X750
Qdot_x = h_g_x.*(T_gas-T_hw);
%This Qdot_x must be matched by the regenerative cooling from gas(rayliegh
%flow)
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


% Resets Matlabs Path preference so it doesn't mess up your matlab
clear CEApath INPPath funcPath

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED
%%
%Chamber Cooling Rayliegh

Lstar = 1.143; %meters, 45 in
chan_ID = 0;
chan_t = 0;
At = (mdot*CStar)/P0; %m^2
contrac = 3; %Empirical chamber contratiuon Ratio for Liquid-Liquid (3 in throat diameter)

[chamber_L, contract_L, nozzle_L, cham_chan_num] = getChamberSize(At, contrac, expansionRatio, Lstar, chan_ID, chan_t);


%Lc = 0.254; %Chamber length based on Huazel empirical correlations

diamt = 2*sqrt(At/pi);
diamt_inches = diamt*39.3701;

Dc = diamt * contrac;
De = diamt * expansionRatio;
num_channels = 20; %iterate this variable?
Chamber_circum = pi*Dc;

A = ((At ./ M_x) .* (((2 + (specificHeatRatio - 1) .* M_x .^ 2) ./ (specificHeatRatio + 1)) .^ ((specificHeatRatio + 1) ./ (2 .* (specificHeatRatio - 1)))));

x = linspace(-(chamber_L + contract_L), nozzle_L, 3200);


% Plot the cross-section
figure()
plot(x, sqrt(A/pi), 'LineWidth', 2);
grid on
yline(0)
xlabel('Radial distance (m)');
ylabel('y (m)');
title('Rocket Nozzle and Chamber Cross-Section');
axis equal

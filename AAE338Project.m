%AAE338 Project

close all;
clear;
clc;

%% Rocket Engine Initialization
pressureExit = 1; %psia
pressureChamber = 400; %psia
OF = 2.35;
mdot_engine = 4; %kg/s propellant mass flow rate (iterate this variable?)
contractionRatio = 3;
Lstar = 1.143; %meters, 45 in
T_hw = 1088; %K -Inconel X750 Wall Temperature (1500 F)
tensile = 5.537*10^8; %Pa - Tensile strength at 1088K (1500 F) 
k = 23; %W/m-K -Nozzle wall made from Inconel X750 @ 1500F
converge_num = 10000; %points in converging section
diverge_num = 10000; %points in diverging section

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
P0 = pressureChamber * 6894.76; %Pa
T_cns = combustionTemperature; %K

R = P0 / (rho0 * T_cns);

Aratio_sub = linspace(contractionRatio, 1, converge_num);
Aratio_sup = linspace(1.001, expansionRatio, diverge_num);
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
% Channel Initial Conditions
helpress_init = 1000000; %Pa
chan_ID = 0.01; %m (1 cm)
saf_fac = 2;
chan_t = saf_fac*(helpress_init*chan_ID/(2*tensile));

At = (mdot_engine * CStar) / P0; %m^2

[chamber_L, contract_L, nozzle_L, cham_chan_num] = getChamberSize(At, contractionRatio, expansionRatio, Lstar, chan_ID, chan_t);

Dt = 2 * sqrt(At/pi);
Dc = sqrt((At * contractionRatio) / pi) * 2;
De = sqrt((At * expansionRatio) / pi) * 2;

A = ((At ./ M_x) .* (((2 + (gma - 1) .* M_x .^ 2) ./ (gma + 1)) .^ ((gma + 1) ./ (2 .* (gma - 1)))));

x1 = linspace(-contract_L, 0, converge_num);
x2 = linspace(0, nozzle_L, diverge_num);
xplot = [x1 x2];
x3 = linspace(-chamber_L, -contract_L, 5);
x4 = sqrt(A(1)/pi) * ones(length(x3));




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

%% Helium Initial Conditions/Loop
%Initial Helium Pressure
Chambertemp = T_cns; %K
heltemp_init = 200; %K
helpress_init = 1000000;%Pa
helmach_init = 0.3;
h_gas = h_g_x(1);
k_wall = k;
wall_thick = 2*chan_t;
Dh = chan_ID;
Cp_init = py.CoolProp.CoolProp.PropsSI("C","T",heltemp_init,"P", helpress_init,"Helium");
Cv_init = py.CoolProp.CoolProp.PropsSI("O","T",heltemp_init,"P", helpress_init,"Helium");
gma_init = Cp_init/Cv_init;
Tstag_hel_init = heltemp_init * (1+((gma_init-1)/2)*helmach_init);
T0stari = heltemp_init * ((1 + gma_init * helmach_init^2)^2 / (2 * (gma_init + 1) * helmach_init^2));

%% Chamber Cooling Loop

% For Loop Conditions
step = 0.01;
Tmax = 1088; % Temperature at which material properties fall apart
Mmax = 1; % Max Mach number allowed in code
steps = floor(2*chamber_L / step); % Number of steps in the foor loop

% Array initialization
qdot_arr = zeros(1, steps);
T_cw_arr = zeros(1, steps);
T_hw_arr = zeros(1, steps);
T_hel_arr = zeros(1, steps);
P_hel_arr = zeros(1, steps);
M_hel_arr = zeros(1, steps);
mdot_arr = zeros(1, steps);

% Informs user Loop is Running
disp('Simulation Running...');

for i = 1:steps
    
    %Setting Helium Step Input Parameters
    if i == 1
        % Sets initial properties of helium
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
    
    %Initialiing Ralyeigh Flow
    Qdot = q_dot * step * (Dh+2*wall_thick);
    rho_i = py.CoolProp.CoolProp.PropsSI("D","T",Ti,"P", Pi,"Helium");
    mdot = rho_i * Vel_i * (pi*(Dh/2)^2);
    
    % Rayleigh Flow Calculations
    [T0e] = getTempStagNew(Qdot, mdot, Mi, Ti, gma_hel, Cp);
    [Me,Te, Pe] = RayleighFlow(Pi, T0e, Mi, gma_hel, T0stari);
    
    % Places helium values into array
    T_hel_arr(i) = Ti;
    P_hel_arr(i) = Pi;
    M_hel_arr(i) = Mi;
    mdot_arr(i) = mdot;

    % Storing Forced Convection Results
    qdot_arr(i) = q_dot;
    T_cw_arr(i) = T_cw;
    T_hw_arr(i) = T_hw;

    % Break Conditions
    Tbreak = T_hw > Tmax;
    MmaxBreak = Me > Mmax;
    
    % Checks to see if the Mach Number has decreased
    if i > 1
        MdecBreak = Me < M_hel_arr(i-1);
    else
        MdecBreak = false;
    end
    
    % Breaks loop if conditions are met
    breakLoop = Tbreak || MmaxBreak || MdecBreak;
    if breakLoop
        % Displays data that could break loop
        disp('Max Value was reached and Loop Was Broken');
        disp('Final Values:');
        fprintf('Hot Wall Temperature:    %0.3f\n', T_hw)
        fprintf('Mach Number:             %0.3f\n', Me);
        if i ~= 1
            fprintf('Previous iteration Mach: %0.3f\n', M_hel_arr(i-1));
        else
            disp('Loop Broke on first iteration.');
        end

        % Breaks Loop
        break
    end
end

% Prints that the simulation is complete
disp('Simulation Complete');

%% Nozzle Cooling Loop

% For Loop Conditions
step_down = chan_ID + 2*chan_t;
Tmax = 1088; % Temperature at which material properties fall apart
Mmax = 1; % Max Mach number allowed in code
steps_down = floor(2*contract_L / step); % Number of steps in the foor loop

% Array initialization
qdot_arr_cha = zeros(1, steps);
T_cw_arr_cha = zeros(1, steps);
T_hw_arr_cha = zeros(1, steps);
T_hel_arr_cha = zeros(1, steps);
P_hel_arr_cha = zeros(1, steps);
M_hel_arr_cha = zeros(1, steps);

% Informs user Loop is Running
disp('Simulation Running...');

for j = 1:steps_down

    [hgas,area,Tgas,tubelen] = nozzleprops(chan_ID, A, h_g_x,...
        chan_t, j,T_gas,converge_num,contract_L,diverge_num, nozzle_L);

    step_around = 0.01;
    steps_around = floor(tubelen/step_around);

    h_gas = hgas;
    Chambertemp = Tgas;

    for i = 1:steps_around
    
        %Setting Helium Step Input Parameters
        if i == 1
            % Sets initial properties of helium
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
    
        %Initialiing Ralyeigh Flow
        Qdot = q_dot * step * (Dh+2*wall_thick);
        rho_i = py.CoolProp.CoolProp.PropsSI("D","T",Ti,"P", Pi,"Helium");
        mdot = rho_i * Vel_i * (pi*(Dh/2)^2);
    
        % Rayleigh Flow Calculations
        [T0e] = getTempStagNew(Qdot, mdot, Mi, Ti, gma_hel, Cp);
        [Me,Te, Pe] = RayleighFlow(Pi, T0e, Mi, gma_hel, T0stari);
    
        % Places helium values into array
        T_hel_arr_cha(i) = Ti;
        P_hel_arr_cha(i) = Pi;
        M_hel_arr_cha(i) = Mi;
    
        % Storing Forced Convection Results
        qdot_arr_cha(i) = q_dot;
        T_cw_arr_cha(i) = T_cw;
        T_hw_arr_cha(i) = T_hw;
    
    end

end

%% Bottom of Script
% Resets Matlabs Path preference so it doesn't mess up your matlab
clear CEApath INPPath funcPath

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED
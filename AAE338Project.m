%% Information

% Class: AAE 338
% Names:
%   - Matthew Currie
%   - Spence Troop
%   - Jonah Fouts
%   - Alex Porter
%   - Preetham Gowni
% Assignment: Final Project

%% Initialization
close all;
clear;
clc;
%lol

%% Graph Settings
chamberGraphs = 1;

%% Value Initialization
% Rocket Engine Initial Conditions
pressureExit = 1; %psia
pressureChamber = 300; %psia
OF = 2.35;
mdot_engine = 0.7; %kg/s propellant mass flow rate (iterate this variable?)
contractionRatio = 3;
Lstar = .4; %meters, 9.84 in
T_hw = 1088; %K -Inconel X750 Wall Temperature (1500 F)
tensile = 5.537*10^8; %Pa - Tensile strength at 1088K (1500 F) 
k_wall = 23; %W/m-K -Nozzle wall made from Inconel X750 @ 1500F
converge_num = 10000; %points in converging section
diverge_num = 10000; %points in diverging section

% Channel Initial Conditions
chan_ID = 0.003; %m (3 mm)
wall_thick = 0.001;
cham_chan_loops = 1; %Number of switch backs chamber has

%Initial Helium Conditons
heltemp_init = 120;
helpress_init = 80e6;
helmach_init = 0.4;
Cp_init = py.CoolProp.CoolProp.PropsSI("C","T",heltemp_init,"P", helpress_init,"Helium");
Cv_init = py.CoolProp.CoolProp.PropsSI("O","T",heltemp_init,"P", helpress_init,"Helium");
gamma_init = Cp_init/Cv_init;

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

[Isp, CStar, expansionRatio, specificHeatRatio, combustionTemperature, Cpcea, ~, ~, rho0,Prcea,visccea,dataeq] = PSP_1DOF_CEA_function_wrapper(pressureChamber,pressureExit, OF, nameString, 0);
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

%%
Pr = (4 * gma) / (9 * gma - 5);
r = Pr ^ (0.33);

T_gas = T_cns .* (1 + (r .* ((gma - 1) ./ 2) .* M_x)) ./ (1 + (((gma - 1) ./ 2) .* M_x));

%% Chamber
At = (mdot_engine * CStar) / P0; %m^2

[chamber_L, contract_L, nozzle_L, cham_chan_num] = getChamberSize(At, contractionRatio, expansionRatio, Lstar, chan_ID, wall_thick);

Dt = 2 * sqrt(At/pi);
Dc = sqrt((At * contractionRatio) / pi) * 2;
De = sqrt((At * expansionRatio) / pi) * 2;
Ac = pi* (Dc/2)^2;

A = ((At ./ M_x) .* (((2 + (gma - 1) .* M_x .^ 2) ./ (gma + 1)) .^ ((gma + 1) ./ (2 .* (gma - 1)))));

x1 = linspace(-contract_L, 0, converge_num);
x2 = linspace(0, nozzle_L, diverge_num);
xplot = [x1 x2];
x3 = linspace(-chamber_L, -contract_L, 5);
x4 = sqrt(A(1)/pi) * ones(length(x3));

hbartz = hgcalc(gma, visccea, Prcea, Cpcea, A./At, pressureChamber * 0.0689476, CStar, Dt);
%h_g_x = (rho_x .* V_x) .^ 0.8; %Convection Coefficent Correlation
%h_g_x = hbartz/h_g_x(1)*h_g_x;

Qdot_x = hbartz .* (T_gas - T_hw);

Thrust = mdot_engine*V_x(end) + pressureExit*A(end);

%% Chamber Cooling Loop

%Initial Conditons for Chamber Loop
h_gas = hbartz(1);
Tstag_hel_init = heltemp_init * (1+((gamma_init-1)/2)*helmach_init);
T0stari = heltemp_init * ((1 + gamma_init * helmach_init^2)^2 / (2 * (gamma_init + 1) * helmach_init^2));


% For Loop Conditions
step = 0.001;
Tmax = 1088; % Temperature at which material properties fall apart
Mmax = 1; % Max Mach number allowed in code
steps = floor(2*cham_chan_loops*chamber_L / step); % Number of steps in the foor loop

% Array initialization
qdot_arr = zeros(1, steps);
T_cw_arr = zeros(1, steps);
T_hw_arr = zeros(1, steps);
T_hel_arr = zeros(1, steps);
P_hel_arr = zeros(1, steps);
M_hel_arr = zeros(1, steps);
mdot_arr = zeros(1, steps);
rho_arr = zeros(1, steps);
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
    [q_dot, T_cw, T_hw] = convergeTemp(T_cns, h_gas, k_wall, wall_thick, Ti, Vel_i, chan_ID, Pi);
    
    %Checking Forced Convection Convergence
    if q_dot == 1
       
        disp('Forced Convection Failed to Converge')
        break
    end
    

    %Initialiing Ralyeigh Flow
    Qdot = q_dot * step * (chan_ID+2*wall_thick);
    rho_i = py.CoolProp.CoolProp.PropsSI("D","T",Ti,"P", Pi,"Helium");
    mdot = rho_i * Vel_i * (pi*(chan_ID/2)^2);
    
    % Rayleigh Flow Calculations
    [T0e] = getTempStagNew(Qdot, mdot, Mi, Ti, gma_hel, Cp);
    [Me,Te, Pe] = RayleighFlow(Pi, T0e, Mi, gma_hel, T0stari);
    
    % Places helium values into array
    T_hel_arr(i) = Ti;
    P_hel_arr(i) = Pi;
    M_hel_arr(i) = Mi;
    mdot_arr(i) = mdot;
    rho_arr(i) = rho_i;
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
step_down = chan_ID + 2*wall_thick;
step_around = 0.001;

Tmax = 1088; % Temperature at which material properties fall apart
Mmax = 1; % Max Mach number allowed in code
steps_down = floor((contract_L + nozzle_L)/ step_down); % Number of steps in the foor loop

max_rad = sqrt(A(end)/pi);
max_tubelen = 2*pi*(max_rad + wall_thick + chan_ID/2);
max_steps_around = floor(max_tubelen/step_around);

% Array initialization
qdot_arr_cha = zeros(steps_down, max_steps_around);
T_cw_arr_cha = zeros(steps_down, max_steps_around);
T_hw_arr_cha = zeros(steps_down, max_steps_around);
T_hel_arr_cha = zeros(steps_down, max_steps_around);
P_hel_arr_cha = zeros(steps_down, max_steps_around);
M_hel_arr_cha = zeros(steps_down, max_steps_around);
rho_arr_cha = zeros(steps_down, max_steps_around);

% Informs user L0oop is Running
disp('Simulation Running...');

for j = 1:steps_down

    [hgas,area,Tgas,tubelen] = nozzleprops(chan_ID, A, hbartz,...
        wall_thick, j,T_gas,converge_num,contract_L,diverge_num, nozzle_L);


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
        [q_dot, T_cw, T_hw] = convergeTemp(Chambertemp, h_gas, k_wall, wall_thick, Ti,Vel_i,chan_ID,Pi);

        %Checking Forced Convection Convergence
        if q_dot == 1

            disp('Forced Convection Failed to Converge')
            break
        end

        %Initialiing Ralyeigh Flow
        Qdot = q_dot * step * (chan_ID+2*wall_thick);
        rho_i = py.CoolProp.CoolProp.PropsSI("D","T",Ti,"P", Pi,"Helium");
        mdot = rho_i * Vel_i * (pi*(chan_ID/2)^2);

        % Rayleigh Flow Calculations
        [T0e] = getTempStagNew(Qdot, mdot, Mi, Ti, gma_hel, Cp);
        [Me,Te, Pe] = RayleighFlow(Pi, T0e, Mi, gma_hel, T0stari);

        % Places helium values into array
        T_hel_arr_cha(j,i) = Ti;
        P_hel_arr_cha(j,i) = Pi;
        M_hel_arr_cha(j,i) = Mi;
        rho_arr_cha(j,i) = rho_i;
        % Storing Forced Convection Results
        qdot_arr_cha(j,i) = q_dot;
        T_cw_arr_cha(j,i) = T_cw;
        T_hw_arr_cha(j,i) = T_hw;

        % Break Conditions
        Tbreak = T_hw > Tmax;
        MmaxBreak = Me > Mmax;

        % Checks to see if the Mach Number has decreased
        if i > 1
            MdecBreak = Me < M_hel_arr_cha(i-1);
        else
            MdecBreak = false;
        end

        % Breaks loop if conditions are met
        breakLoop = Tbreak || MmaxBreak || MdecBreak;
        if breakLoop
            % Displays data that could break loop
            disp('Max Value was reached and Loop Was Broken');
            disp('Final Values:');
            fprintf('Hot Wall Temperature: %0.3f\n', T_hw)
            fprintf('Mach Number:          %0.3f\n', Me);
            if i ~= 1
                fprintf('Previous Mach:    %0.3f\n', M_hel_arr_cha(i-1));
            else
                disp('Loop Broke on first iteration.');
            end

            % Breaks Loop
            break
        end

    end

    percentDone = j / steps_down * 100;
    fprintf('%0.2f Percent Complete\n', percentDone);
end
%% Thrust Calc and Plotting
mdot_hel_total = (ceil(cham_chan_num)/(cham_chan_loops*2) + steps_down)*mdot;

fprintf("Engine Thrust:         %0.2f N \n", Thrust)
fprintf("Helium Mass Flow Rate: %0.2f kg/s \n", mdot_hel_total)

disp('Simulation Complete');

M_hel_arr_cha(M_hel_arr_cha==0) = NaN;
xscale = [0:step_around:(size(M_hel_arr_cha,2)-1)*step_around];
yscale = [0:step_down:(size(M_hel_arr_cha,1)-1)*step_down];
figure()
s = pcolor(xscale,yscale,M_hel_arr_cha);
colorbar
s.EdgeColor = 'none';
title('Helium Mach Number Distribution in Nozzle Loop')
xlabel('Circumfrential Distance [m]')
ylabel('Distance from End of Chamber [m]')
a=colorbar;
ylabel(a,'Mach Number')

P_hel_arr_cha(P_hel_arr_cha==0) = NaN;
figure()
s = pcolor(xscale,yscale,P_hel_arr_cha./1000);
s.EdgeColor = 'none';
title('Helium Pressure Distribution in Nozzle Loop')
a=colorbar;
ylabel(a,'Pressure [kPa]')
xlabel('Circumfrential Distance [m]')
ylabel('Distance from End of Chamber [m]')

qdot_arr_cha(qdot_arr_cha==0) = NaN;
figure()
s = pcolor(xscale,yscale,qdot_arr_cha);
s.EdgeColor = 'none';
a=colorbar;
title('Steady State Heat Flux Distribution in Nozzle Loop')
ylabel(a,'Heat Flux [W/m^2]')
xlabel('Circumfrential Distance [m]')
ylabel('Distance from End of Chamber [m]')

rho_arr_cha(rho_arr_cha==0) = NaN;
figure()
s = pcolor(xscale,yscale,rho_arr_cha);
s.EdgeColor = 'none';
a=colorbar;
title('Helium Density Distribution in Nozzle Loop')
ylabel(a,'Desnity [kg/m^3]')
xlabel('Circumfrential Distance [m]')
ylabel('Distance from End of Chamber [m]')

T_hw_arr_nozz = T_hw_arr_cha;
T_hw_arr_nozz(T_hw_arr_nozz==0) = NaN;
figure()
s = pcolor(xscale,yscale,T_hw_arr_nozz);
s.EdgeColor = 'none';
a=colorbar;
title('Nozzle Wall Temperature Distribution in Nozzle Loop')
ylabel(a,'Temperature [K]')
xlabel('Circumfrential Distance [m]')
ylabel('Distance from End of Chamber [m]')


% Plots Graphs
if chamberGraphs
    coolinggrapher(M_hel_arr, qdot_arr, T_hw_arr, T_cw_arr, T_hel_arr,...
        P_hel_arr, Aratio, xplot, A, T_gas, hbartz, Qdot_x, M_x, rho_x,...
        T_x, V_x, x3, x4, T_hw_arr_cha, M_hel_arr_cha, P_hel_arr_cha, qdot_arr_cha, rho_arr_cha, T_hw_arr_nozz, rho_arr, hbartz)

end


%% Bottom of Script
% Resets Matlabs Path preference so it doesn't mess up your matlab
clear CEApath INPPath funcPath

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED




%% Referenced Functions
% This example includes the external functions (No CEA program files)
%
% <include>funcs/RayleighFlow.m</include>
%
% <include>funcs/nozzleprops.m</include>
%
% <include>funcs/getTempStagNew.m</include>
%
% <include>funcs/getRe.m</include>
%
% <include>funcs/getHc.m</include>
%
% <include>funcs/getChamberSize.m</include>
%
% <include>funcs/coolinggrapher.m</include>
%
% <include>funcs/convergeTemp.m</include>
%
% <include>funcs/PSP_1DOF_CEA_function_wrapper.m</include>
%
% The file content above is properly syntax highlighted.






%% LOOP

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

    step_around = 0.01;
    steps_around = floor();
    
    h_gas = h_g_x();
    Chambertemp = t_gas();

    for i = 1:steps_around
    
        %Setting Helium Step Input Parameters
        if i == 1
            % Sets initial properties of helium
            Ti = heltemp_init;
            Pi = helpress_init;
            Mi = helmach_init;
            T0i = Tstag_hel_init;
    
            % Places intial values into array
    %         T_hel_arr(1) = Ti;
    %         P_hel_array(1) = Pi;
    %         M_hel_arr(1) = Mi;
    
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
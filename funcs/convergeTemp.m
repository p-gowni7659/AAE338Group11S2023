function [q_dot, T_cw, T_hw] = convergeTemp(T_gas, h_gas, k, wall_thick,...
    T_i,Vel,Dh,P_i)

% Grabs helium convection heat transfer coefficient
h_hel = getHc(Vel,Dh,T_i,P_i);

% Solves system of equations
num = (h_gas / h_hel) * T_gas + (wall_thick / k) * h_gas * T_gas + T_i;
den = (h_gas / h_hel) + 1 + (wall_thick / k) * h_gas;

% Derives other values
T_hw = num / den;
q_dot = h_gas * (T_gas - T_hw);
T_cw = -(q_dot * wall_thick / k) + T_hw;

end

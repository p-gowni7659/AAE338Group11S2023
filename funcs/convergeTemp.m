function [q_dot, T_cw, T_hw] = convergeTemp(T_gas, h_gas, k, wall_thick, T_i,mdot,Dh,P_i)

T_hw_guess = T_gas - 100;


tolerance = 10;
step = 1;

trials = 100000;
i = 1;

T_hw = T_hw_guess;

while i < trials

    q_dot_hc = h_gas*(T_gas - T_hw);

    T_cw = -(q_dot_hc * wall_thick/ k) + T_hw;

    h_hel = getHc(mdot,Dh,T_i,P_i,T_cw);

    q_dot_cc = h_hel*(T_cw - T_i);

    if abs(q_dot_hc - q_dot_cc) <= tolerance

        q_dot = q_dot_cc;
        i = trials;

    elseif q_dot_hc > q_dot_cc

        T_hw = T_hw + step;
        q_dot = 1;

    else

        T_hw = T_hw - step;
        q_dot = 1;

    end

    i = i + 1;

end
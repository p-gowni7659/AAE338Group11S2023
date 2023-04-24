function [q_dot, T_cw, T_hw] = convergeTemp(T_gas, h_gas, k, wall_thick, T_i)

T_hw_guess = T_gas;

tolerance = 0.1;

trials = 1000;
i = 1;

T_hw = T_hw_guess;

while i <= trials

    q_dot_hc = h_gas*(T_gas - T_hw);

    T_cw = -(q_dot_hc * wall_thick/ k) + T_hw;

    q_dot_cc = h_cg*(T_cw - T_i);

    if abs(q_dot_hc - q_dot_cc) <= tolerance

        q_dot = q_dot_cc;

    else

        T_hw = (T_hw + T_cw)/2;

    end

    i = i + 1;

end
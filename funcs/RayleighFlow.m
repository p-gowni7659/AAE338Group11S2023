function [Me,Te, Pe] = RayleighFlow(Pi, T0e, Mi, gamma, T0star)
    % Calculates New Mach Number using Stagnation Temperature of the 
    % exit flowing gas

    % Calculates T0* of the flowing gas and ratio
    T_ratio = T0e / T0star;

    % Solves for Mach Number
    frac1 = (-gamma * T_ratio) / (gamma^2*T_ratio - gamma^2 + 1);
    
    frac2num = sqrt(-4*gamma^2*T_ratio - 8*gamma*T_ratio - 4*T_ratio +...
        4*gamma^2 + 8*gamma + 4);
    frac2den = 2 * (gamma^2 * T_ratio - gamma^2 + 1);
    frac2 = frac2num / frac2den;

    frac3 = (gamma + 1) / (gamma^2*T_ratio - gamma^2 + 1);

    Me = abs(sqrt(frac1 - frac2 + frac3));

     % Calculates Static Presure at Exit
    Pe = Pi * ((gamma * Mi^2 + 1) / (gamma * Me^2 + 1));

    % Calculates Static Temperature at Exit
    Te = T0e / (1 + ((gamma - 1) / 2) * Me^2);

end

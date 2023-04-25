function [Me,Te, Pe] = RayleighFlow(Pi, Ti, T0e, Mi, gamma)
    % Calculates New Mach Number using Stagnation Temperature of the 
    % exit flowing gas

    % Calculates T0* of the flowing gas and ratio
    T0star = Ti * ((1 + gamma * Mi^2)^2 / (2 * (gamma + 1) * Mi^2));
    T_ratio = T0e / T0star;
    % T0i = Ti * (1+ ((gamma-1)/2)*Mi^2);
    %disp(T0star)
    % Creates Equation to solve for
    syms Mo
    Mnum = 2 * (gamma + 1) * Mo^2;
    Mden = (1 + gamma * Mo^2)^2;
    Mmult = 1 + ((gamma - 1) / 2) * Mo^2;
    M_func = T_ratio == (Mnum / Mden) * Mmult;

    % Solves for Mach Number and Exit Static Temperature
    MeArray = abs(double(solve(M_func, Mo)));
    disp(MeArray)
    Me = MeArray(MeArray > 0 & MeArray <= 1);

    % If critical condition is met, Mach 1 number is taken
    if length(Me) > 1
        Me = Me(1);
    end

    % Calculates Static Presure at Exit
    Pe = Pi * ((gamma * Mi^2 + 1) / (gamma * Me^2 + 1));

    % Calculates Static Temperature at Exit
    Te = T0e / (1 + ((gamma - 1) / 2) * Me^2);

end

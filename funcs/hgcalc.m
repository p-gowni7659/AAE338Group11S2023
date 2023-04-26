function [hg] = hgcalc(gma, visc, Pr, Cp, A_At, Pcns, cstar, Dt)
Dt_R = 1.084145663;
sig = 1.4;
hg = (41.8565/(Dt^0.2))  * (((visc^0.2)*Cp)/(Pr^0.6)) * (((Pcns*9.81)/cstar)^0.8) * (Dt_R^0.1) * ((1/A_At)^0.9) * sig;
end
function [Te0] = getTempStagNew(Qdot, mdot, Mi, Ti, Pi, gamma, Cp)
% Calculates T0e after the application of the heat
% Grabs Cp value and Calculates Exit Temperature
Ti0 = Ti * (1 + ((gamma - 1) / 2) * Mi^2);
Te0 = (Qdot / (mdot * Cp)) + Ti0;
end

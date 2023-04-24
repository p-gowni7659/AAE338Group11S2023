function [Re] = getRe(mdot, hydraulic_D, dynamicViscosity)
% Returns the reynolds number

Re = (4 * mdot) / (pi * hydraulic_D * dynamicViscosity);

end
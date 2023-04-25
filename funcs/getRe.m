function [Re] = getRe(Vel, dynvisc, rho, Dh)
% Returns the reynolds number

Re = (rho*Vel*Dh)/dynvisc;

end
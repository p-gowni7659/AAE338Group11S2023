function [Hc] = getHc(vel,Dh,Ti,Pi)
dy_vsic_bulkT = py.CoolProp.CoolProp.PropsSI('V','T',Ti,'P',Pi,'Helium');
rho = py.CoolProp.CoolProp.PropsSI('D','T',Ti,'P',Pi,'Helium');
Pr = py.CoolProp.CoolProp.PropsSI('Prandtl','T',Ti,'P',Pi,'Helium');
k_bulkT = py.CoolProp.CoolProp.PropsSI('L','T',Ti,'P',Pi,'Helium');
Re = (rho*vel*Dh)/dy_vsic_bulkT;

Nu = 0.023*(Re^0.8)*(Pr^(0.4));
Hc = (k_bulkT/Dh)*Nu;
end
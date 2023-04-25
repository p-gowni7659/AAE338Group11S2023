function [Hc] = getHc(mdot,Dh,Ti,Pi,Tcw)
dy_vsic_bulkT = py.CoolProp.CoolProp.PropsSI('V','T',Ti,'P',Pi,'Helium');
dy_vsic_surfT = py.CoolProp.CoolProp.PropsSI('V','T',Tcw,'P',Pi,'Helium');
Pr = py.CoolProp.CoolProp.PropsSI('Prandtl','T',Ti,'P',Pi,'Helium');
k_bulkT = py.CoolProp.CoolProp.PropsSI('L','T',Ti,'P',Pi,'Helium')*1000;
Re = getRe(mdot,Dh,dy_vsic_bulkT);
Nu = 0.027*(Re^0.8)*(Pr^(0.4));
Hc = (k_bulkT/Dh)*Nu;
end
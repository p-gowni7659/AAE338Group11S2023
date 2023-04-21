%AAE338 Project
%Combustion Gases - RP1/Lox, mixture ratio 2.35, Molar Mass 10.20583
%kg/mol
gma = 1.222;
P0 = 6894.76;
T_cns = 3588.889;
rho0 = P0/(R*T_cns);
R = 8.314/10.20583;
%As an example, used F1 Engine Epxanison Ratio of 16:1. Guessed Chamber
%Area Ratio of 3:1 since it will change.
Aratio_sub = linspace(3,1,200);
Aratio_sup = linspace(1.01,16,3000);
Aratio = [Aratio_sub,Aratio_sup];
M_x = [];
rho_x = [];
T_x = [];
for subratio = Aratio_sub
    [mach_sub,T_sub,P_sub,rho_sub,area_sub] = flowisentropic(gma,subratio,'sub');
    M_x = [M_x mach_sub];
    rho_x = [rho_x rho_sub];
    T_x = [T_x T_sub];
end
for supratio = Aratio_sup
    [mach_sup,T_sup,P_sup,rho_sup,area_sup] = flowisentropic(gma,supratio,'sup');
    M_x = [M_x mach_sup];
    rho_x = [rho_x rho_sup];
    T_x = [T_x T_sup];
end
rho_x = rho_x .* rho0;
T_x = T_x .* T_cns;
V_x = M_x .* sqrt(gma.*R.*T_x);
plot(Aratio, M_x)
grid on;
title('Nozzle Mach Number Distribution')
xlabel("Area Ratio [A/At]")
ylabel('Mach Number')
figure()
plot(Aratio, rho_x)
grid on;
title('Nozzle Density Distribution')
xlabel("Area Ratio [A/At]")
ylabel('Density [kg/m^3]')
figure()
plot(Aratio, T_x)
grid on;
title('Nozzle Temperature Distribution [Isentropic]')
xlabel("Area Ratio [A/At]")
ylabel('Temperature [K]')
figure()
plot(Aratio, V_x)
grid on;
title('Nozzle Velocity Distribution')
xlabel("Area Ratio [A/At]")
ylabel('Velocity [m/s]')
%%
Pr = (4*gma)/(9*gma - 5);
r = Pr^(0.33);

T_gas = T_cns .* (1+(r.*((gma-1)./2).*M_x))./(1+(((gma-1)./2).*M_x));

h_g_x = (rho_x .* V_x).^0.8;
T_hw = 1473.15; %K -Inconcel X750 Wall Temperature
k = 25; %W/m-K -Nozzle wall made from Inconel X750
Qdot_x = h_g_x.*(T_gas-T_hw);
%This Qdot_x must be matched by the regenerative cooling from gas(rayliegh
%flow)
figure()
plot(Aratio, T_gas)
grid on;
title('Nozzle Temperature Gas Distribution [Correlation]')
xlabel("Area Ratio [A/At]")
ylabel('Gas Temperature [K]')
figure()
plot(Aratio, Qdot_x)
grid on;
title('Nozzle Qdot Distribution')
xlabel("Area Ratio [A/At]")
ylabel('Qdot')
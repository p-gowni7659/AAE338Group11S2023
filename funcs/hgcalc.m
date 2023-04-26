function [hgarrc] = hgcalc(gma, visc, Pr, Cp, A_At, Pcns, cstar, Dt)
Dt_R = 1.084145663;
sig = 0.3;
sigarr = ones(1,10000) .* sig;
sigarr = [sigarr linspace(0.3,0.25,10000)];

hgarr = (41.8565./(Dt.^0.2))  .* (((visc.^0.2).*Cp)./(Pr.^0.6)) .* (((Pcns.*9.81)./cstar).^0.8) .* (Dt_R.^0.1) .* ((1./A_At).^0.9) .* sigarr;

Rd =  0.001;

hgarrc = (1./(1./hgarr) + Rd);
end
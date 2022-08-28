function res = crf1ex
% Membrane separation process for binary feed using cross-flow model
clear all;
crfdata;
Pa = Pm(1); Pb = Pm(2); D1 = ((1-alpa)*r + alpa)/2;
F1 = -((1-alpa)*r - 1)/2; E1 = alpa/2 - D1.*F1;
R1 = 1./(2*D1-1); S1 = (alpa.*(D1-1) + F1)./((2*D1 - 1).*(alpa/2 - F1));
T1 = 1./(1 - D1 - E1./F1);
i0 = xr./(1-xr); i2 = xf./(1-xf);
ur = -D1.*i0 + sqrt((D1.^2).*i0.^2 + 2*E1.*i0 + F1.^2);
uf = -D1.*i2 + sqrt((D1.^2).*i2.^2 + 2*E1.*i2 + F1.^2);
theta = 1 - ((1-xf)./(1-xr)).*((uf-E1./D1)./(ur-E1./D1)).^R1...
.*((uf-alpa+F1)./(ur-alpa+F1)).^S1 .*((uf-F1)./(ur-F1)).^T1;
yp = (xf - (1-theta)*xr)/theta;
Ami = quad(@(x) mArea(x,D1,E1,F1,R1,S1,T1,alpa,r,xf,uf), i0, i2);
Am = Ami*qf*t/(ph*Pb); rc = theta*yp./xf; % recovery ratio
% Calculated variables: yp(permeate mole fraction), Am(area)
% theta(stage-cut) or xr(reject composition)
% Results: results = [yp, xr, Am, theta, rc]
res.yp = yp; res.theta = theta; res.Am = Am; res.rc = rc; res.xr = xr;
end

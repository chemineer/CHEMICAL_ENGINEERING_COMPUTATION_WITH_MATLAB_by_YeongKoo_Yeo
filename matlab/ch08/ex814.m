% heatH2Ocal.m: calculation of H2O heater
clear all; heatH2Odat; b0 = 0.0;
b = fsolve(@heatH2O,b0);
Ai = Qr/(Ui*F12*dTlm); % inside heat transfer area
Np = ceil(Ai/(At*Nt)); % number of tube passes
Acf = Ds*cl*b/Pt; % cross-sectional area between baffles and shell axis
Go = mh/Acf; % shell side mass velocity
Nreo = De*Go/(muh*3600/1488); % shell-side Reynolds number
Npro = Cph*(muh*3600/1488)/kh; %  shell-side Prandtl number
Nuo = 0.36*Nreo^0.55*Npro^(1/3)*(muh/muwh)^0.14; % shell-side Nusselt number
ho = Nuo*kh/De; % shell-side heat transfer coefficient
Nrei = Di*rhoc*ui/(muc/1488);  % tube-side Reynolds number
Npri = Cpc*(muc*3600/1488)/kc; % tube-side Prandtl number
Nui = 0.027*Nrei^(0.8)*Npri^(1/3)*(muc/muwc)^0.14; % tube-side Nusselt number
hi = Nui*kc/Di; % tube-side heat transfer coefficient
Ui = 1/((Di/Do)/ho + (Di*dx)/(Dlm*kw) + 1/hi + Rfi + (Di/Do)*Rfo);
Uo = Ui*Di/Do;
xt = Pt/Do; xl = xt; Ks = 1.1*L/b; Nr = ceil(Ds/Pt/2);
fp = (0.044 + 0.08*xl/((xt - 1)^(0.43+1.13/xl)))*((Do*Go)/muh)^(-0.15);
dPo = 2*Ks*Nr*fp*(Go/3600)^2/(gc*rhoh*144); % shell-side pressure drop (psi)
fD = 1/(1.82*log10(Nrei) - 1.64)^2; Gi = rhoc*ui;
dPi = 0.6*Np*fD*Gi^2*L/(gc*rhoc*Di*144); % tube-side pressure drop (psi)
Ai = Qr/(Ui*F12*dTlm); % inside heat transfer area
Ao = Qr/(Uo*F12*dTlm); % outside heat transfer area
fprintf('Number of tube passes = %g\n',Np);
fprintf('Baffle spacing = %g in\n',b*12);
fprintf('Shell-side overall heat transfer coefficient = %g Btu/(ft^2*h*F)\n', Uo);
fprintf('Tube-side overall heat transfer coefficient = %g Btu/(ft^2*h*F)\n',Ui);
fprintf('Tube inside heat transfer area = %g ft^2\n',Ai);
fprintf('Tube outside heat transfer area = %g ft^2\n',Ao);
fprintf('Shell-side pressure drop = %g psi\n',dPo);
fprintf('Tube-side pressure drop = %g psi\n',dPi);

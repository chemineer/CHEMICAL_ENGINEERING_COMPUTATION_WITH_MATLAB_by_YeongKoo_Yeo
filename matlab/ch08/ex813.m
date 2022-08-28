% coolcal.m: calculation of heat exchanger
coolhxdat; b0 = 0.0; % guess initial b
b = fsolve(@coolhxf,b0);
Acf = Ds*cl*b/Pt; % cross-sectional area between baffles and shell axis
Go = mc/Acf; % shell side mass velocity
Nreo = De*Go/(muc*3600/1488); % shell-side Reynolds number
Npro = Cpc*(muc*3600/1488)/kc; %  shell-side Prandtl number
Nuo = 0.36*Nreo^0.55*Npro^(1/3)*(muc/muwc)^0.14; % shell-side Nusselt number
ho = Nuo*kc/De; % shell-side heat transfer coefficient
Nrei = Di*rhoh*ui/(muh/1488);  % tube-side Reynolds number
Npri = Cph*(muh*3600/1488)/kh; % tube-side Prandtl number
Nui = 0.027*Nrei^(0.8)*Npri^(1/3)*(muh/muwh)^0.14; % tube-side Nusselt number
hi = Nui*kh/Di; % tube-side heat transfer coefficient
Ui = 1/((Di/Do)/ho + (Di*dx)/(Dlm*kw) + 1/hi + Rfi + (Di/Do)*Rfo); % new overall heat transfer coefficients
Uo = Ui*Di/Do; xt = Pt/Do; xl = xt; Ks = 1.1*L/b; Nr = ceil(Ds/Pt/2);
fp = (0.044 + 0.08*xl/((xt - 1)^(0.43+1.13/xl)))*((Do*Go)/muc)^(-0.15);
dPo = 2*Ks*Nr*fp*(Go/3600)^2/(gc*rhoc*144); % shell-side pressure drop (psi)
fD = 1/(1.82*log10(Nrei) - 1.64)^2; Gi = rhoh*ui;
Ai = Q/(Ui*F12*dTlm); Ao = Q/(Uo*F12*dTlm);
Np = ceil(Ai/(At*Nt)); % number of tube passes
dPi = 0.6*Np*fD*Gi^2*L/(gc*rhoh*Di*144); % tube-side pressure drop (psi)
fprintf('Baffle spacing = %g in\n',b*12);
fprintf('Shell-side overall heat transfer coefficient = %g Btu/ (ft^2*h*F)\n',Uo);
fprintf('Tube-side overall heat transfer coefficient = %g Btu/ (ft^2*h*F)\n',Ui);
fprintf('Tube inside heat transfer area = %g ft^2\n',Ai);
fprintf('Tube outside heat transfer area = %g ft^2\n',Ao);
fprintf('Shell-side pressure drop = %g psi\n',dPo);
fprintf('Tube-side pressure drop = %g psi\n',dPi);

function fun = coolhxf(b)
% coolhxn.m: design of condenser
coolhxdat; % retireve data
Ai = Q/(Ui*F12*dTlm); % inside heat transfer area
Np = ceil(Ai/(At*Nt)); % number of tube passes
Ai = Np*Nt*pi*Di*L; % adjust Ai based on Np
De = (4/(pi*Do))*(Pt^2 - pi*Do^2/4); % hydraulic effective diameter
Acf = Ds*cl*b/Pt; % cross-sectional area between baffles and shell axis
Go = mc/Acf; % shell-side mass velocity
Nreo = De*Go/(muc*3600/1488); % shell-side Reynolds number
Npro = Cpc*(muc*3600/1488)/kc; % shell-side Prandtl number
Nuo = 0.36*Nreo^0.55*Npro^(1/3)*(muc/muwc)^0.14; % shell-side Nusselt number
ho = Nuo*kc/De; % shell-side heat transfer coefficient
Nrei = Di*rhoh*ui/(muh/1488); % tube-side Reynolds number
Npri = Cph*(muh*3600/1488)/kh; % tube-side Prandtl number
Nui = 0.027*Nrei^(0.8)*Npri^(1/3)*(muh/muwh)^0.14; % tube-side Nusselt number
hi = Nui*kh/Di; % tube-side heat transfer coefficient
Ui = 1/((Di/Do)/ho + (Di*dx)/(Dlm*kw) + 1/hi + Rfi + (Di/Do)*Rfo);
Qd = Ai*Ui*F12*dTlm; fun = Q - Qd;
end

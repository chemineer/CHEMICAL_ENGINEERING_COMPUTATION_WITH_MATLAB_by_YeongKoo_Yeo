function fun = heatH2O(b)
% heatH2O.m: design of H2O heater
heatH2Odat;
Ai = Qr/(Ui*F12*dTlm); % inside heat transfer area
Np = ceil(Ai/(At*Nt)); % number of tube passes
Ai = Np*Nt*pi*Di*L; % adjust Ai based on Np
Ui = Qr/(Ai*F12*dTlm); % adjust Ui
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
Qd = Ai*Ui*F12*dTlm;
fun = Qr - Qd;
end

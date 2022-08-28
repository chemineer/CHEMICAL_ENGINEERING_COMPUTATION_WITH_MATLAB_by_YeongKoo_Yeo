% heatH2Odat.m
% Physical properties (hot(shell): demineralized water, cold(tube): raw water)
Cph = 1.0; Cpc = 1.01; rhoh = 62.4; rhoc = 62.4; kh = 0.36; kc = 0.363;
muh = 0.81; muc = 0.92; muwc = muc; muwh = muh;
% Operating conditions
mh = 1.5e5; ui = 5; T1 = 95; T2 = 85; t1 = 75; t2 = 80;
% Shell and tube
L = 10; Do = 0.0625; Di = 0.05167; Ds = 1.77083; Pt = 1/12;
cl = 0.02083; Rfi = 0.001; Rfo = 0; kw = 30;
% Parameters and basic properties
gc = 32.2; Qr = mh*Cph*(T1 - T2); % heat load (Btu/h)
mc = Qr/(Cpc*(t2 - t1)); % cold stream rate (lb/h)
dTlm = ((T1-t2) - (T2-t1))/log((T1-t2)/(T2-t1)); % log-mean temp.
R = (T1-T2)/(t2-t1); S = (t2-t1)/(T1-t1);
F12den = (R-1)*log((2 - S*(R+1-sqrt(R^2+1)))/(2 - S*(R+1+sqrt(R^2+1))));
F12 = sqrt(R^2+1)*log((1-S)/(1-R*S))/F12den; % correction factor
Aci = mc/(rhoc*ui*3600); % tube-side total cross-sectional area (ft^2/pass)
Nt = ceil(4*Aci/(pi*Di^2)); % number of tubes per pass
At = pi*Di*L; % heat transfer area per tube
dx = (Do - Di)/2; % tube thickness
Dlm = (Do - Di)/log(Do/Di); % log-mean diameter
De = (4/(pi*Do))*(Pt^2 - pi*Do^2/4); % hydraulic effective diameter
Ui = 400; % Assume Ui

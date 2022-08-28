% coolhxdat.m
% Physical properties
Cph = 0.92; Cpc = 1; rhoh = 59.76; rhoc = 59.87; kh = 0.3; kc = 0.36;
kw = 30; muh = 0.75; muc = 0.77; muwc = muc; muwh = muh;
% Operating conditions
mh = 6.2e4; ui = 5; T1 = 150; T2 = 120; t1 = 75; t2 = 100;
% Shell and tube
L = 15; Do = 0.0625; Di = 0.04017; Ds = 1.4375; Pt = 1/12; cl = 0.02083; Rfi = 0.004;
Rfo = 0;
% Parameters and basic properties
gc = 32.2; dTlm = ((T1-t2) - (T2-t1))/log((T1-t2)/(T2-t1)); % log-mean temp.
R = (T1-T2)/(t2-t1); S = (t2-t1)/(T1-t1);
F12den = (R-1)*log((2 - S*(R+1-sqrt(R^2+1)))/(2 - S*(R+1+sqrt(R^2+1))));
F12 = sqrt(R^2+1)*log((1-S)/(1-R*S))/F12den; % correction factor
Q = mh*Cph*(T1 - T2); mc = Q/(Cpc*(t2 - t1)); % Heat load and cold stream rate
Aci = mh/(rhoh*ui*3600); % tube-side total cross-sectional area
Nt = ceil(4*Aci/(pi*Di^2)); % number of tubes per pass
At = pi*Di*L; % heat transfer area per tube
dx = (Do - Di)/2; % tube thickness
Dlm = (Do - Di)/log(Do/Di); % log-mean diameter
De = (4/(pi*Do))*(Pt^2 - pi*Do^2/4); % hydraulic effective diameter
Ui = 150; % Assume Ui

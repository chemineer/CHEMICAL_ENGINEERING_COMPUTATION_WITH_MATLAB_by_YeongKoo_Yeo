function DPs = dpshell(D,Ds,L,Lbc,Lbin,Lbout,Lc,Dotl,Dsb,Pt,Nt,Nss,w,rho,Visc,Phi,Layout)
% Calculate shell-side pressure drop
% input:
%  D: tube outside diameter (mm)     Ds: shell inside diameter(mm)
%  L: tube length (m)                Lbc: central baffle spacing (mm)
%  Lbin: inlet baffle spacing (mm)   Lbout: outlet baffle spacing (mm)
%  Lc: baffle cut (mm)               Dotl: shell outer tube limit (mm)
%  Dsb: shell-baffle clearance (mm)  Pt: tube pitch (mm)
%  Nt: total number of tubes in the bundle  Nss: number pairs of sealing strips
%  w- shell-side flow rate (kg/s)    rho: fluid density (kg/m^3)
%  Visc: shell-side fluid viscosity (Ns/m^2)    Phi: viscosity correction factor
%  Layout: tube layout (1:triangular, 2:in-line square, 3:rotated square)
% output:
%  DPs: shell-side pressure drop
% Correlational coefficients for tube arrangement
b1 = [0.372 0.486 4.570 45.100 48.000; 0.391 0.0815 6.090 32.100 35.000;
0.303 0.333 3.500 26.200 32.000];
b2 = -[0.123 0.152 0.476 0.973 1.000; 0.148 -0.022 0.602 0.963 1.000; 0.126 0.136 0.476 0.913 1.000];
b3 = [7.00 6.30 6.59]; b4 = [0.500 0.378 0.520];
switch Layout
case 1, Pp = 0.866*Pt; Pn = Pt/2; Pd = Pt;
case 2, Pp = Pt; Pn = Pt; Pd = Pn;
otherwise, Pp = 0.7071*Pt; Pn = Pp; Pd = Pn;
end
% Friction factor for ideal tube-bank
Nc = Ds*(1 - 2*Lc./Ds)/Pp; Ptd = Pt/D; Sm = Lbc*(Ds - Dotl + (Pt-D)*(Dotl-D)/Pd);
Nres = D*w*1e3/(Visc*Sm); % Shell-side Reynolds number
if Nres >= 1e4, J = 1; c1 = 3.7; elseif Nres >= 1e3, J = 2; c1 = 3.7;
elseif Nres >= 100, J = 3; c1 = 3.7; elseif Nres >= 10, J = 4; c1 = 4.5; else J = 5; c1 = 4.5; end
b = b3(Layout)/(1 + 0.14*Nres^b4(Layout)); Fj = b1(Layout,J) * (1.33/Ptd)^b * Nres^b2(Layout,J);
Dpbi = 2*Fj*Nc*(w/(Sm*1e-6))^2 / (rho*Phi);
% Pressure drop for ideal window section
bdm = (Ds - 2*Lc)/Dotl; bd1 = acos(bdm); bd2 = 1 - 2*Lc/Ds;
Ftc = (pi + 2*bdm*sin(bd1) - 2*bd1)/pi;
Sw = (Ds^2/4)*(acos(bd2)-bd2*sqrt(1-bd2^2)) - (Nt/8)*(1-Ftc)*pi*D^2;
Dw = 4*Sw/(1.5708*Nt*(1-Ftc)*D + 2*Ds*bd2); Ncw = 0.81*Lc/Pp; Gw = 1e6*w/ sqrt(Sm*Sw);
if Nres >= 100, Dpw = 0.5*Gw^2 * (2+0.6*Ncw)/rho;
else, Dpw = (2.6e4*Visc*(Ncw/(Pt-D) + Lbc/Dw^2) + Gw)*Gw/rho; end
% Correction factor for baffle leakage effects
Stb = 0.6223*D*(1 + Ftc)*Nt; Ssb = Ds*Dsb*0.5*(pi - acos(bd2));
R1 = (Stb + Ssb)/Sm; R2 = Ssb/(Ssb + Stb); Pv = -0.15*(1+R2) + 0.8;
Rl = exp(-1.33*(1 + R2)*R1^Pv);
% Correction factor for bundle bypassing
Fbp = (Ds - Dotl)*Lbc/Sm; Nsc = Nss/Nc;
if Nsc >= 0.5, Rb = 1;
else
if Nss == 0, c2 = 0; else c2 = (2*Nsc)^0.3333; end
Rb = exp(-c1*c2*Fbp);
end
% Correction factor due to unequal baffle spacing at inlet and outlet
Nb = 1e3*L/Lbc + 1;
if Nres >= 100, An = 0.2; else An = 1; end
Rs = (Lbin/Lbc)^(2-An) + (Lbout/Lbc)^(2-An);
% Shell-side pressure drop
DPs = Rb*Dpbi*((Nb-1)*Rl + (1+Ncw/Nc)*Rs) + Nb*Rl*Dpw;
end

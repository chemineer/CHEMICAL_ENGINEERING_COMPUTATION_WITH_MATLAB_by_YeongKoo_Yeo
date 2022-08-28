function Hshell = htcshell(D,Ds,L,Lbc,Lbin,Lbout,Lc,Dotl,Dsb,Pt,Nt,Nss,w,Visc,Cp,Xk,Phi,Layout)
% Calculate shell-side heat transfer coefficient
% input:
%  D: outside diameter of tube (mm)     Ds: inside diameter of shell (mm)
%  L: tube length (m)           Lbc: central baffle spacing (mm)
%  Lbin: inlet baffle spacing (mm)      Lbout: outlet baffle spacing (mm)
%  Lc: baffle cut (mm)          Dotl: shell outer tube limit(mm)
%  Dsb: shell-baffle clearance (mm)     Pt: tube pitch (mm)
%  Nt: total number of tubes in the bundle  Nss: number of pairs of sealing strips
%  w- flow rate in shell (kg/s)     Visc: shell-side fluid viscosity (Ns/m^2)
%  Cf: heat capacity of shell-side fluid (J/kg/K)   Xk: heat conductivity of shell-side fluid (W/m/K)
%  Phi: viscosity correction factor
%  Layout: tube layout (1:triangular, 2:in-line square, 3:rotated square)
% output:
%   Hshell: shell-side heat transfer coefficient
%   Correlational coefficients for tube arrangement
a1 = [0.321 0.321 0.593 1.360 1.400; 0.370 0.107 0.408 0.900 0.970; 0.370 0.370 0.730 0.498 1.550];
a2 = -[0.388 0.388 0.477 0.657 0.667; 0.395 0.266 0.460 0.631 0.667; 0.396 0.396 0.500 0.656 0.667];
a3 = [1.450 1.187 1.930]; a4 = [0.519 0.370 0.500];
switch Layout
case 1, Pp = 0.866*Pt; Pn = Pt/2; Pd = Pt;
case 2, Pp = Pt; Pn = Pt; Pd = Pn;
otherwise, Pp = 0.7071*Pt; Pn = Pp; Pd = Pn;
end
Nc = Ds*(1 - 2*Lc./Ds)/Pp; Sm = Lbc*(Ds - Dotl + (Pt-D)*(Dotl-D)/Pd);
Nres = D*w*1e3/(Visc*Sm); % ShellÃø Reynolds ¼ö
Pr = Cp*Visc/Xk; Ptd = Pt/D;
% Heat transfer coefficients for ideal tube bank
% (Nres: Reynolds number of shell-side fluid)
if Nres >= 1e4, J = 1; c1 = 1.25;
elseif Nres >= 1e3, J = 2; c1 = 1.25;
elseif Nres >= 100, J = 3; c1 = 1.25;
elseif Nres >= 10, J = 4; c1 = 1.35;
else J = 5; c1 = 1.35; end
a = a3(Layout)/(1 + 0.14*Nres^a4(Layout)); Hj = a1(Layout,J) * (1.33/Ptd)^a * Nres^a2(Layout,J);
Hsi = Hj*Cp*Phi*w*Pr^(-0.6667) / (Sm*1e-6);
% Correction factor for baffle configuration effects
adm = (Ds - 2*Lc)/Dotl; adm1 = acos(adm); Ftc = (pi + 2*adm*sin(adm1) - 2*adm1)/pi;
Phic = Ftc + 0.54*(1 - Ftc)^0.345;
% Correction factor for baffle leakage effects
Stb = 0.6223*D*(1 + Ftc)*Nt; Ssb = Ds*Dsb*0.5*(pi - acos(1 - 2*Lc/Ds));
R1 = (Stb + Ssb)/Sm; R2 = Ssb/(Ssb + Stb); Phil = 0.44*(1-R2) + (1-0.44*(1-R2))*exp(-2.2*R1);
% Correction factor for bundle bypassing
Fbp = (Ds - Dotl)*Lbc/Sm; Nsc = Nss/Nc;
if Nsc >= 0.5, Phib = 1;
else, if Nss == 0, c2 = 0; else c2 = (2*Nsc)^0.3333; end
Phib = exp(-c1*c2*Fbp);
end
Nb = 1e3*L/Lbc + 1;
% Correction factor for adverse temperature gradient buildup
if Nres >= 100, Phir = 1;
else
Ncw = 0.8*Lc/Pp; Phs = 1.51/((Nc+Ncw)*(Nb+1))^0.18;
if Nres <= 20, Phir = Phs;
elseif Nres <= 100, Phir = Phs - (1-Phs)*(0.25-0.0125*Nres); end
if Phir <= 0.4, Phir = Phs; end
end
% Correction factor due to unequal baffle spacing at inlet and outlet
if Nres >= 100, An = 0.6; else An = 0.333; end
Phis = (Nb-1 + (Lbin/Lbc)^(1-An) + (Lbout/Lbc)^(1-An))/(Nb-1 + (Lbin+Lbout)/Lbc);
% Shell-side heat transfer coefficient
Hshell = Hsi*Phic*Phil*Phib*Phir*Phis;
end

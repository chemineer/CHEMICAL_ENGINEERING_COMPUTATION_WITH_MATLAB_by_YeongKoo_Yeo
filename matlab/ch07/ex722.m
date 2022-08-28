% mshortcut.m: multicomponent distillation by shortcut method
% 1: n-hexane(C6H14), 2: n-heptane(C7H16), 3: n-octane(C8H18)
% LK: n-hexane(1), HK: n-heptane(2), HNK: n-octane(3)
clear all;
% Data
xf = [0.33, 0.37, 0.3]; xd = [0.99, 0.01, 0]; xb(1) = 0.01; % mole fraction
F = 100; Tf = 105; P = 1.2; % F: mol/h, T: deg.C, P: atm
q = 0.4; % feed stream is 60% vapor
% Antoine parameters
A = [6.87024, 6.89385, 6.90940]; B = [1168.720, 1264.370, 1349.820];
C = [224.210, 216.636, 209.385];
% Determine product flow rates and composition
% z(1)=D, z(2)=xb(2), z(3)=xb(3)
mf = @(z) [F*xf(1)-z(1)*xd(1)-(F-z(1))*xb(1); F*xf(2)-z(1)*xd(2)-(F-z(1))*z(2);
F*xf(3)-z(1)*xd(3)-(F-z(1))*z(3)];
z0 = [5, 0.5, 0.5]; z = fsolve(mf,z0); Dist = z(1); xb(2) = z(2); xb(3) = z(3);
Bott = F - Dist; % bottom product rate
% LK and HK compositions
xlkf = xf(1); xhkf = xf(2); xlkd = xd(1); xhkd = xd(2); xlkb = xb(1); xhkb = xb(2);
% Determine boiling points of top stream (D) and bottom stream (B)
Df = @(Td) (760*P)*sum(xd./10.^(A - B./(Td + C))) - 1;
Bf = @(Tb) sum(xb.*10.^(A - B./(Tb + C)))/(760*P) - 1;
Td0 = 100; Tb0 = 100;
Td = fzero(Df,Td0); % boiling point of distillate (D)
Tb = fzero(Bf,Tb0); % boiling point of bottom product (B)
% Find distribution coefficients Ki = yi/xi = Pivp/P
Kf = 10.^(A - B./(Tf + C))/(760*P); % feed stream
Kd = 10.^(A - B./(Td + C))/(760*P); % distillate
Kb = 10.^(A - B./(Tb + C))/(760*P); % bottom product stream
% Relative volatility (alpa_LK/HK = alpa12) and average alpa
alpaf12 = Kf(1)/Kf(2); % feed stream
alpad12 = Kd(1)/Kd(2); % distillate
alpab12 = Kb(1)/Kb(2); % bottom product stream
avgalpa = (alpaf12*alpad12*alpab12)^(1/3); % average relative volatility
% (1) Minimum number of trays
Nmin = log((xlkd/xhkd)*(xhkb/xlkb))/log(avgalpa); % Fenske eqn.
% Minimum reflux ratio: Underwood eqn.
alpaf = Kf/Kf(2); % relative volatility of feed stream
alpad = Kd/Kd(2); % relative volatility of distillate
funf = @(theta) 1 - q - sum(alpaf.*xf./(alpaf - theta)); theta0 = 1.5;
theta = fzero(funf,theta0);
Rmin = sum(alpad.*xd./(alpad - theta)) - 1; % minimum reflux ratio
Ract = 1.5*Rmin; % R = 1.5*Rmin
% (2) ideal number of trays if R=1.5Rmin: Gilliland correlation
X = (Ract - Rmin)/(Ract + 1); Y = 0.75*(1 - X^0.5668); % Eduljee correlation
N = ceil((Nmin + Y)/(1 - Y)); % number of trays (integer)
% (3) feed tray location
Nrs = (xhkf*xlkb^2*Bott/(xlkf*xhkd^2*Dist))^0.206; % Kirkbride eqn.
Ns = ceil(N/(Nrs + 1)); % number of stripping section trays
Nr = N - Ns; % number of rectifying section trays
Ftray = Nr; % feed tray location from top
% Results
fprintf('Dew point(D) = %g (deg.C), bubble point(B) = %g (deg.C)\n', Td, Tb);
fprintf('(1) the minimum number of trays = %g\n', Nmin);
fprintf('(2) the number of ideal trays (R=1.5*Rmin) = %g\n', N);
fprintf('(3) feed tray location = %g\n', Ftray);

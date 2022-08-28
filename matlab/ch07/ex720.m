% fugdist.m
% 1: ethane, 2: propane(LK), 3: i-butane(HK), 4: n-butane
xF = [0.15 0.18 0.18 0.49]; % feed composition (mole fraction)
A = [3.93835 4.53678 4.3281 4.35576]; B = [659.739 1149.36 1132.108 1175.581];
C = [-16.719 24.906 0.918 -2.071]; P = 7;
% Compositions of D and B (D: distillate, B: bottom)
Br(1) = 0; Br(3) = 0.95*xF(3); Br(4) = xF(4); Br(2) = 0.001*(Br(3) + Br(4))/0.999;
Dr(1) = xF(1); Dr(2) = xF(2) - Br(2); Dr(3) = 0.05*xF(3); Dr(4) = 0;
xD = Dr/sum(Dr); xB = Br/sum(Br);
% Dew point(D) and bubble point(B)
Td0 = 270; Tb0 = 300;
fD = @(Td) sum(P*xD./(10.^(A - B./(Td + C)))) - 1;
fB = @(Tb) sum(xB.*10.^(A - B./(Tb + C))/P) - 1;
Td = fzero(fD, Td0); Tb = fzero(fB, Tb0);
% Fenske equation
kD = (10.^(A - B./(Td + C)))/P; kB = (10.^(A - B./(Tb + C)))/P;
num = log((xD(2)/xD(3))*(xB(3)/xB(2))); den = log(sqrt((kD(2)/kD(3)*(kB(2)/ kB(3)))));
Nm = num/den;
% Temperature of feed stream
Tf0 = 300; fF = @(Tf) sum(xF.*10.^(A - B./(Tf + C))/P) - 1;
Tf = fzero(fF, Tf0); kF = (10.^(A - B./(Tf + C)))/P;
% Underwood equation
th0 = 1.5;
fTh = @(th) sum(xF.*10.^(A - B./(Tf + C))/P./(10.^(A - B./(Tf + C))/P - th*10^(A(3) - B(3)/(Tf+C(3)))/P));
theta = fzero(fTh, th0); alpD = kD/kD(3); Rm = sum(alpD.*xD./(alpD - theta)) - 1;
% Gilliland equation
R = 1.5*Rm; X = (1.5 - 1)*Rm/(R + 1); Y = 0.75*(1 - X^0.5658); N = (Y+Nm)/(1-Y);
% Kirkbride equation
c1 = sum(Br)/sum(Dr); c2 = xF(3)/xF(2); c3 = (xB(2)/xD(3))^2;
h = (c1*c2*c3)^0.206; p = N/(1+h); m = p*h;
% Print results
fprintf('Bubble point(B) = %8.4f, Dew point(D) = %8.4f', Tb, Td);
fprintf('\nFeed temp.(F) = %8.4f', Tf);
fprintf('\ntheta = %8.4f, min. number of stages = %7.4f', theta, Nm);
fprintf('\nMin. reflux ratio = %7.4f, actual number of stages = %7.4f',N, Rm);
fprintf('\nStages above the feed stage = %7.4f', m);
fprintf('\nStages below the feed stage = %7.4f\n', p);

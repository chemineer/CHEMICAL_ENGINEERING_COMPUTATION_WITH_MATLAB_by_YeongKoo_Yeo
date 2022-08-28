% hxnst.m: Shell-and-tube heat exchanger
clear all;
hxndat; % input data
% Problem type
if ptype > 0
if ptype <= 3, Nvar = 1;     %

calculate tube-side inlet temp.(Ti1)
elseif ptype <= 6, Nvar = 2;   %   calculate tube-side outlet temp.(Ti2)
else Nvar = 3; end       %calculate tube-side flow rate(Wi)

end
Rw = 1e-3*(Do - Di)/(2*Xkw);  % Xkw: inverse of tube-wall heat conductivity(1/hw)
A = pi*1e-3*Do*(L - 2e-3*Ls)*Nt;
if ptype == 0
Cpim = hxcp(cpreft,Trt,(Ti1+Ti2)/2);  % tube side
Cpsm = hxcp(cprefs,Trs,(Ts1+Ts2)/2);  % shell side
Qi = Wi*Cpim*(Ti2 - Ti1);  Qs = Ws*Cpsm*(Ts2 - Ts1);
Qr = abs(Qi/Qs);
if (abs(1-Qr) > 0.1)
display('Tube-side and shell-side heat duty differ by more than 10%.');
display('Tube-side heat duty is used as total duty.');
end
Q = Qi;
end
% Iterations until differences in successive values of key variables converge to zero.
varC = 10; iter = 0;
while varC > 1e-3
if ptype > 0
% Calculation of unknown variables using energy balances
Cpim = hxcp(cpreft,Trt,(Ti1+Ti2)/2);  % tube-side
Cpsm = hxcp(cprefs,Trs,(Ts1+Ts2)/2);  % shell-side
Qs = Ws*Cpsm*(Ts2 - Ts1); Qi = -Qs;
iter1 = 0;
switch Nvar
case 1  % calculates Ti1(tube-side inlet temperature)
Ti1new =  Ti2 - Qi/(Wi*Cpim); crT = 10;
while crT >= 1e-3
Ti1 = Ti1new; Cpim = hxcp(cpreft,Trt,(Ti1+Ti2)/2);
Ti1new = Ti2 - Qi/(Wi*Cpim); crT = abs((Ti1new - Ti1)/Ti1new);
iter1 = iter1 + 1;
end
case 2  %  calculates Ti2(tube-side outlet temperature)
Ti2new = Ti1 + Qi/(Wi*Cpim); crT = 10;
while crT >= 1e-3
Ti2 = Ti2new; Cpim = hxcp(cpreft,Trt,(Ti1+Ti2)/2);
Ti2new = Ti1 + Qi/(Wi*Cpim); crT = abs((Ti2new - Ti2)/Ti2new);
iter1 = iter1 + 1;
end
case 3  % calculates Wi(tube-side flow rate)
Wi = abs(Qs/(Cpim*(Ti2 - Ti1)));
end
end
Tib = (Ti1 + Ti2)/2;
% Shell-side heat transfer coefficient
Tsb = (Ts1 + Ts2)/2;
mus = hxvis(murefs,Trs,Tsb,fsS); % shell-side viscosity
rhos = hxrho(rhorefs,Trs,Tsb,fsS); % shell-side density
Cps = hxcp(cprefs,Trs,Tsb); % shell-side heat capacity
xks = hxthc(xkrefs,Trs,Tsb); % shell-side heat conductivity
Tw = (Tsb + Tib)/2; Twnew = Tw; crT = 10;
while crT >= 1e-3
if fsS == 1 % liquid
Phis = (mus/hxvis(murefs,Trs,Tw,1))^0.14;
elseif fsS == 2  % gas
Phis = (Tsb/Tw)^0.25;
end
Hs = htcshell(Do,Ds,L-2e-3*Ls,Lbc,Lbin,Lbout,Lc,Dotl,Dsb,Pt,Nt,Nss,Ws,mus,Cps,xks,Phis,Layout);
% tube-side heat conductivity
mut = hxvis(mureft,Trt,Tib,fsT); % tube-side viscosity
rhot = hxrho(rhoreft,Trt,Tib,fsT); % tube-side density
Cpt = hxcp(cpreft,Trt,Tib); % tube-side heat capacity
xkt = hxthc(xkreft,Trt,Tib); % tube-side heat conductivity
Ui = 4e6*Wi*Npass/(pi*rhot*Di^2*Nt); Rei = 1e-3*Di*Ui*rhot/mut;
Pri = Cpt*mut/xkt;
if fsT == 1, Phit = (mut/hxvis(mureft,Trt,Tw,1))^0.14;
elseif fsT == 2, Phit = (Tw/Tib)^0.25; end
Ht = htctube(Rei,Pri,Di,L,xkt,Phit);
% Estimate tube wall temperature
Tw = Tib + Hs/(Hs + Ht) * (Tsb - Tib); crT = abs((Tw - Twnew)/Tw); Twnew = Tw;
end  % end while

% Calculate heat duty from heat transfer equations
U = 1/(Do/(Di*Ht) + 1/Hs + Rw + Rds + Rdt); Dt1 = Ts1 - Ti2; Dt2 = Ts2 - Ti1;
if (Dt1 <= 0 || Dt2 <= 0) break; end
Delt = (Dt1 - Dt2)/log(Dt1/Dt2); Ft = 1;
if Npass > 1
R = (Ts1 - Ts2)/(Ti2 - Ti1); P = (Ti2 - Ti1)/(Ts1 - Ti1); tm = sqrt(R^2 + 1);
Ft = tm*log((1-P)/(1-R*P))/((R-1)*log((2-P*(R+1-tm))/(2-P*(R+1+tm))));
end
Deltm = Delt*Ft;
if ptype == 0, Areq = Qi/(U*Deltm); Da = (A - Areq)/A*100; break;
end
Q = U*A*Deltm;
% Calculate assumed unknown variable
Sgn = 1; if Qs < 0, Sgn = -1; end
switch ptype
case {1,4,7}
Ts1new = Ts2 - Sgn*Q/(Ws*Cps); varC = abs((Ts1 - Ts1new)/Ts1new); Ts1 = Ts1new;
case {2,5,8}
Ts2new = Ts1 + Sgn*Q/(Ws*Cps); varC = abs((Ts2 - Ts2new)/Ts2new); Ts2 = Ts2new;
case {3,6}
Wsnew = abs(Q/((Ts2 - Ts1)*Cps)); varC = abs((Ws - Wsnew)/Wsnew); Ws = Wsnew;
end
iter = iter + 1;
end
% Pressure drop
DPs = dpshell(Do,Ds,L-2e-3*Ls,Lbc,Lbin,Lbout,Lc,Dotl,Dsb,Pt,Nt,Nss,Ws,rhos,mus,Phis,Layout);
DPt = dptube(Di,L,rf,Ui,rhot,mut,Phit,Npass);
% Print results
fprintf('Overall heat transfer coefficient: U = %g(W/m^2/K)\n', U);
fprintf('Heat transfer coefficient: tube-side = %g(W/m^2/K), shell-side = %g(W/m^2/K)\n',Ht,Hs);
fprintf('Heat duty: Q = %g(W)\n', Q);
fprintf('Pressure drop: tube-side = %g(Pa), shell-side = %g(Pa)\n', DPt,DPs);
fprintf('Tube-side: Ti1 = %g(K), Ti2 = %g(K), flow rate = %g(kg/sec)\n',Ti1,Ti2,Wi);
fprintf('Shell-side: Ts1 = %g(K), Ts2 = %g(K), flow rate = %g(kg/sec)\n',Ts1,Ts2,Ws);

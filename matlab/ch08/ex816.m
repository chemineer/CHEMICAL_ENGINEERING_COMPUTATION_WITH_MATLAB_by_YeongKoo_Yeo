% htcoef.m : heat transfer coefficients (double-pipe finned-tube heat exchanger)
clear all;
% Data
Fft = 0.002; Kt = 0.366; mut = 0.72; Cpt = 1; Jt = 9.75; % Tube
Ffs = 0.002; Ks = 0.074; mus = 2.45; Cps = 0.518; % Shell
Dis = 3.068/12; Dos = 3.5/12; % shell size
Dit = 1.61/12; Dot = 1.9/12; % tube size
N = 24; Hf = 0.5/12; w = 0.035/12; % fin size
Gt = 26740; Gs = 18000; % flow rates
% Tube-side :
Deq = Dit; Pr = Cpt*(2.4191*mut)/Kt; % 1 cP = 2.4191 lb/(ft-hr)
Re = Gt*Deq/(2.4191*mut);
if (Re < 10000), hit = Jt* Pr^0.333 * Kt/Deq;
else, hit = 0.023 * Re^0.8 * Pr^0.333 * Kt / Deq; end
hift = 1 / (1/hit + Fft);
% Shell-side (bare tube):
Deq = Dis - Dot; Pr = Cps*(2.4191*mus)/Ks; % 1 cP = 2.4191 lb/(ft-hr)
Re = Gs*Deq/(2.4191*mus);
if (Re < 10000), his = Jt* Pr^0.333 * Ks/Deq;
else, his = 0.023 * Re^0.8 * Pr^0.333 * Ks / Deq; end
hifs = 1 / (1/his + Ffs);
% Shell-side (finned-tube):
Af = 2*Hf*N; Cs = pi*(Dis^2 - Dot^2)/4; Nf = Cs - w*Af/2;
Deq = 4*Nf/(pi*(Dis+Dot) - N*w + Af); Ao = pi*Dot + Af; X = Hf*sqrt(hifs/(6*Kt*w));
e = (exp(X) - exp(-X))/(exp(X) + exp(-X))/X; ep = e*Af/Ao + 1 - Af/Ao; hifd = hifs*ep;
Ar = pi*Dit/(pi*Dot + Af); Uo = 1/((Dot-Dit)/Kt + hift*Ar + hifd);
% Output
fprintf('Shell-side heat transfer coefficient: %10.7f\n', hifd);
fprintf('Tube-side heat transfer coefficient: %10.7f\n', hift);
fprintf('Fin efficiency: %10.7f\n', ep);
fprintf('Overall heat transfer coefficient: %10.7f\n', Uo);

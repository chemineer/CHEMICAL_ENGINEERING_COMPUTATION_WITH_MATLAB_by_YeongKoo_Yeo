% optrf.m: optimal reflux ratio in a binary distillation column
clear all;
% Data
opdat.we = 0.8; opdat.eta = 0.75; opdat.S = 12000; opdat.K = 0.05;
opdat.rhoG = 2; opdat.rhoL = 850; opdat.rhos = 8000;
opdat.lamb = 800; opdat.lambs = 1800; opdat.Css = 0.05; opdat.Cst = 10;
opdat.F = 100; opdat.T = 70; opdat.P = 760; opdat.alpa = 2.3;
opdat.xB = 0.05; opdat.xD = 0.85; opdat.xF = 0.4;
% Assign data
alpa = opdat.alpa; xD = opdat.xD; xF = opdat.xF;
% Find optimal reflux ratio rfopt
rf0 = 1.5; % initial guess
rfopt = fminsearch(@rfobj,rf0,[],opdat);
minrf = (xD/xF - alpa*(1-xD)/(1-xF))/(alpa-1); % minimum rf by Fenske eqn.
rfpr = 1.2*minrf; % reflux ratio by the rule of thumb (=1.2*minrf)
fprintf('Optimal reflux ratio = %g\n', rfopt);
fprintf('Reflux ratio by the rule of thumb = %g\n', rfpr);

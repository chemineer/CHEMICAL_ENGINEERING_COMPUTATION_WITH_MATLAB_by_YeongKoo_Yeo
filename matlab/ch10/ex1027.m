% optmm.m: determination of optimal parameters for Michaelis-Menten model
clear all;
S = [0.74, 1.19, 1.46, 1.62, 1.63, 5.27, 5.87, 6.26, 8.32, 10.10, 11.10,...
20.10, 21.73, 25.10, 27.78, 35.70]; % substrate concentration
r = [0.10, 0.21, 0.21, 0.14, 0.34, 0.49, 0.36, 0.48, 0.82, 1.26, 0.55,...
3.34, 2.49, 2.01, 1.80, 1.68]; % reaction rate
rmax0 = 1.7; km0 = 11; % initial guesses
J = @(x) sum((r - x(1)*S./(x(2) + S)).^2); % x(1) = rmax, x(2) = km
x = fminsearch(J,[rmax0, km0]); rmax = x(1); km = x(2); % optimization
fprintf('Maximum reaction rate (rm) = %g\n',rmax);
fprintf('Michaelis constant (km) = %g\n',km);
% Compare data and model
Si = linspace(floor(min(S)),ceil(max(S)),1000);
ri = rmax*Si./(km + Si); % model
plot(Si,ri,S,r,'o'), xlabel('S(ng/ml)'), ylabel('r(nmol/min)')
legend('Michaelis-Menten model','Data','location','best')

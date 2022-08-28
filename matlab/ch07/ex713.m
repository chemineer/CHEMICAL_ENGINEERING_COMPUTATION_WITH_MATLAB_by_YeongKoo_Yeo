% trayabs.m
clear all;
x = [0.0 0.033 0.072 0.117 0.171]; % equilibrium data (x)
ye = [0.0 0.0396 0.0829 0.1127 0.136]; % equilibrium data (y*)
V = 85; % solute-free gas flow rate (kmol/hr)
x0 = 0.002; yp = 0.01; % mole fractions at the top of the column

% (1) Curve fitting on the equilibrium data (order <=4)
P3 = polyfit(x,ye,3);  % fit the data by 3rd-order polynomial
P4 = polyfit(x,ye,4);  % fit the data by 4th-order polynomial
y3 = polyval(P3,x); y4 = polyval(P4,x);
rmse3 = sqrt(sum((ye-y3).^2)); rmse4 = sqrt(sum((ye-y4).^2)); % find RMSE
n = 3; P = P3; if rmse3 > rmse4, n = 4; P = P4; end % choose order of the fitting curve
fprintf('Order of the fitting polynomial = %g\n', n);
xi = 0:0.001:0.18; yi = polyval(P,xi); % generate data for plot
subplot(2,2,1), plot(xi,yi,x,ye,'o'), xlabel('x'), ylabel('y^*')
title('Equilibrium curve'), legend('Fitting curve','Data','location','best')

% (2) plots of operating curves for L = 170, 150 and 130 110 kmol/hr
L = [170 150 130]; r = L/V; ypr = yp/(1-yp); x0r = x0/(1-x0);
y = [];
for k = 1:length(L), w = r(k)*xi./(1-xi) + ypr - r(k)*x0r; z = w./(1 + w); y = [y z']; end
subplot(2,2,2), plot(xi,yi,xi,y(:,1),xi,y(:,2),'--',xi,y(:,3),':')
xlabel('x'), ylabel('y'), title('Operating curves')
legend('Equil. curve','L=170','L=150','L=130','location','best')

% (3) find the point at which the operating curve becomes tangent to the equilibrium curve
% use bisection method
La = 80; Lb = 120; Lm = (La + Lb)/2; crit = abs(La - Lb);
r = La/V; w = r*xi./(1-xi) + ypr - r*x0r; opa = w./(1 + w);
r = Lb/V; w = r*xi./(1-xi) + ypr - r*x0r; opb = w./(1 + w);
r = Lm/V; w = r*xi./(1-xi) + ypr - r*x0r; opm = w./(1 + w);
Da = min(opa - yi); Db = min(opb - yi); Dm = min(opm - yi);
while crit > 1e-6
if (Da*Dm) < 0, Lb = Lm; Db = Dm; Lm = (La + Lb)/2;
else, La = Lm; Da = Dm; Lm = (La + Lb)/2; end
r = Lm/V; w = r*xi./(1-xi) + ypr - r*x0r; opm = w./(1 + w);
Dm = min(opm - yi); crit = abs(La - Lb);
end
Lmin = Lm; fprintf('Minimum liquid flow rate = %g kmol/h\n', Lmin);

% (4) plot of the operating line when L = Lmin
rm = Lmin/V; w = rm*xi./(1-xi) + ypr - rm*x0r; z = w./(1 + w);
subplot(2,2,3),, plot(xi,yi,xi,z,'--')
xlabel('x'), ylabel('y'), title('Operating curve at L=L_{min}')
legend('Equil. curve','L=L_{min}','location','best')
% determine x at which the operating curve becomes tangent to the
% equilibrium curve
f = @(x) (rm*x./(1-x) + ypr - rm*x0r)/(1+rm*x./(1-x) + ypr - rm*x0r) -...
(P(1)*x.^4 + P(2)*x.^3 + P(3)*x.^2 + P(4)*x + P(5));
x0 = 0.03; xmin = fsolve(f,x0);
fprintf('The operating curve becomes tangent to the equilibrium curve at x = %g.\n', xmin);

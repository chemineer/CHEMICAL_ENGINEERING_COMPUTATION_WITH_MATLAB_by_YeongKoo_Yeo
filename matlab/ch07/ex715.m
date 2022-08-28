% absdiam.m: diameter of a SO2 absorber
% Data
rhog = 1.17; rhol = 1000; Fp = 130; phi = 1; mul = 0.8; Lm = 2450; G = 103; g = 9.82; f = 0.7;
% Find diameter
L = 1.5*Lm; X = (L/G)*sqrt(rhog/rhol); z = -1.668 - 1.085*log(X) - 0.297*(log(X))^2; Y = 10^z;
Gf = sqrt(rhog*rhol*g*Y/(phi*Fp*mul^0.2)); % superficial velocity
Gr = 0.7*Gf; A = (G/60)/Gr; Dt = sqrt(4*A/pi);
fprintf('Column diameter = %g m\n', Dt)

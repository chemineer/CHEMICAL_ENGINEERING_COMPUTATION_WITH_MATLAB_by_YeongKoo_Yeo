V = 1; F = 25; Caf = 10; k1 = 50; k2 = 100; k3 = 10; % data
f = @(Cs) [-k1*Cs(1) - k3*Cs(1)^2 + F*(Caf - Cs(1))/V; k1*Cs(1) - k2*Cs(2) - F*Cs(2)/V];
Cs0 = [5 5]; Cs = fsolve(f,Cs0) 
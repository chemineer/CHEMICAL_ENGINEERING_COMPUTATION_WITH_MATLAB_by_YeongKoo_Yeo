function dz = LTmodel(t,z,ht)
% z(1) = h1, z(2) = h2, z(3) = T1, z(4) = T2
% Retrieve data
A1 = ht.A1; A2 = ht.A2; F0 = ht.F0; T0 = ht.T0; H = ht.H;
c1 = ht.c1; c2 = ht.c2; rCp = ht.rCp; Q1 = ht.Q1; Q2 = ht.Q2;
h1 = z(1); h2 = z(2); T1 = z(3); T2 = z(4);
% Flow rates
F1b = c1*sqrt(h1 - h2); F2 = c2*sqrt(h2); dh = max(H,h1) - max(H,h2); F1t = c1*sqrt(dh);
% Heat input
qi1 = Q1/(rCp*A1*h1); qi2 = Q2/(rCp*A2*h2);
% Differential equations
dz = [(F0 - F1t - F1b)/A1; (F1t + F1b - F2)/A2;
(F0/A1/h1)*(T0 - T1) + qi1; (F1t + F1b)*(T1 - T2) + qi2];
end

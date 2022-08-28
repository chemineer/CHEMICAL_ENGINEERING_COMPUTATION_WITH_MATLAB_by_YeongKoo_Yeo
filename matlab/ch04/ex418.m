% btraoult.m: bubble point temperature by Raoult's law 
% mixture: 32% n-hexane(1), 31% n-heptane(2), 25% n-octane(3), 12% n-nonane(4) 
A = [4.00139, 4.02023, 4.05075, 4.07356]; % Antoine constants 
B = [1170.875, 1263.909, 1356.360, 1438.03]; 
C = [224.317, 216.432, 209.635, 202.694]; 
x = [0.32 0.31 0.25 0.12]; % Compositions 
P = 1.5; % total pressure (bar) 
f = @(T) sum(x.*10.^(A-B./(C+T-273.15)))/P - 1; % define the nonlinear equation (T: K) 
T0 = 400; T = fsolve(f,T0); T = T - 273.15 % T0: initial guess 
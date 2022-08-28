% pfraxd.m: PFR with axial diffusion using the method of line (MoL) 
clear all; 
n = 50; pf.n = n; pf.Pe = 1; pf.Da = 2; % Data and parameters  
h = 1/n; Z0 = ones(n,1); tspan = [0 1]; % initialize 
[t,C] = ode45(@pfrdiff,tspan,Z0,[],pf); % solve the set of ODEs 
Cs = [1, C(end,:)]; % steady-state by MoL 
x = [h:h:1]; mesh(x,t,C); xlabel('\xi'), ylabel('\tau'), zlabel('\phi (\tau,\xi)') 
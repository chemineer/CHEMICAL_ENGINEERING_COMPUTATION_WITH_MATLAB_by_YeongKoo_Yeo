% pfrmol.m: calculation of PFR using the method of lines (MoL) 
clear all; 
n  =  20; pf.k  =  0.18; pf.v  =  0.5; pf.C0  =  1; pf.L  =  0.5; % Data and parameters 
pf.n  =  n; h  =  pf.L/n; C0  =  pf.C0; 
Z0  =  ones(n,1)*C0; tspan  =  [0 10]; % initialize 
[t,C]  =  ode45(@pfrconc,tspan,Z0,[],pf); % solve the set of ODEs 
Cs  =  [C0, C(end,:)]; % steady-state by MoL 
x  =  [0:h:pf.L]; Cm  =  C0*exp(-pf.k*x/pf.v); % steady-state by exact solution 
% Plot results 
subplot(1,2,1), plot(x,Cs,x,Cm,'--'), xlabel('x(m)'), ylabel('C(mol/ liter)') 
legend('St-st by MoL','St-st by exact solution') 
subplot(1,2,2), plot(t,C(:,end)), xlabel('t(min)'), ylabel('C(mol/liter)')
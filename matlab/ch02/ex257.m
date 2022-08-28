% daerxn.m: solution of DAE (differential algebraic equations) system 
clear all; 
a  =  200; k  =  0.02; kg  =  0.01; u  =  1; Ca0  =  1; L  =  1; % data 
Cas0  =  kg*Ca0/(k+kg); M  =  [1 0;0 0]; opt  =  odeset('Mass',M); 
zspan  =  [0 L]; y0  =  [Ca0; Cas0]; % initial conditions 
[z,y]  =  ode15s(@hrdae,zspan,y0,opt,a,u,k,kg); 
plot(z,y(:,1),z,y(:,2),'--'), xlabel('Location(z)'), ylabel('Concentration') 
legend('C_A','C_{As}') 
function dy  =  hrdae(z,y,a,u,k,kg) 
% Differential algebraic equations - heterogeneous reactor (1st-order reaction) 
Ca  =  y(1); Cas  =  y(2); 
dy  =  [-(kg*a/u)*(Ca - Cas); % differential equation 
kg*(Ca - Cas) - k*Cas]; % algebraic equation 
end 
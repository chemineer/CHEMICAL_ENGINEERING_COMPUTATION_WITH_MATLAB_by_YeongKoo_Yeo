% qdP.m: plot of volumetric flow rate vs. dP 
clear all; 
dP = -18*1.01325e5; dz = 100; rho = 1e3; mu = 1e-3; L = 1500; D = 0.154; rh = 4.57e-5; % data 
x0 = [5 0.001]; % initial guess (x(1)=v, x(2)=f) 
x = fsolve(@qdpfun,x0,[],dP,dz,rho,mu,L,D,rh); 
v = x(1); f = x(2); Q = x(1)*pi*D^2/4; % volumetric flow rate 
dPf = 2*f*rho*L*v^2/D; % pressure drop due to friction loss 
fprintf('Volumetric flow rate of water = %g m^3/sec\n', Q); 
fprintf('Pressure drop due to friction loss = %g kPa\n', dPf/1000); 

function fun = qdpfun(x,dP,dz,rho,mu,L,D,rh) 
v = x(1); f = x(2); g = 9.8; Nre = D*v*rho/mu; % Reynolds number 
fun = [-v^2/2 + g*dz + dP/rho + 2*f*L*v^2/D; 1/sqrt(f) + 1.7372*log(rh/3.7/D + 1.255/Nre/sqrt(f))]; 
end 
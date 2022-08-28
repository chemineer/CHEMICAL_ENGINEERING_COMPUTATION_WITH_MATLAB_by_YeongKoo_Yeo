function [avgv] = vtwall(rho,mu,delta) 
% Calculates velocity distribution and average velocity for vertical laminar flow of a falling film 
% input 
% rho: density (kg/m^3), mu: viscosity (kg/m/s) 
% delta: film thickness (m)  
% output 
% avgv: average velocity (m/s) 
g = 9.8; h = delta/100; x = [0:h:delta]; Rg = rho*g*delta^2/2/mu; v = Rg.*(1 - (x/delta).^2); avgv = Rg*2/3; 
plot(x,v), grid, ylabel('v_z(x)'), xlabel('x'), axis([0 delta 0 max(v)]) 
end 
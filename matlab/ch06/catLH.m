function z = catLH(y,De,n,h,k,Kr,Ca0) 
% Diffusion and reaction in a catalyst pellet by Langmuir-Hinshelwood kinetics 
% y: dimensionless concentration, De: effective diffusivity 
m = n+1;    
for i = 1:m, x(i) = h*(i-1); end % radial distance 
for i = 1:m    
rxn = k*y(i)/sqrt(1 + Kr*y(i)^2);    
if i == 1, z(i) = 2*De*(y(i+1) - y(i))/h^2 - rxn;    
elseif i == m, z(i) = y(i) - Ca0; % i = n+1    
else % i = 2,3,...,n        
z(i) = De*(y(i+1)-2*y(i)+y(i-1))/h^2 + De*(y(i+1)-y(i-1))/(h*x(i)) - rxn;    
end 
end 
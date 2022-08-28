% denlqbz.m: density of liquid benzene 
Tc = 288.93; Pc = 49.24; w = 0.212;  % Tc: deg.C, Pc: bar 
Mw = 78; % molecular weight (g/mol) 
T = 0:70; n = length(T); 
for k = 1:n          
denBz(k) = gunyam(Tc,Pc,w,Mw,T(k)); % density: kg/m^3 
end 
plot(T,denBz), grid, xlabel('T(deg.C)'), ylabel('Density(kg/m^3)')
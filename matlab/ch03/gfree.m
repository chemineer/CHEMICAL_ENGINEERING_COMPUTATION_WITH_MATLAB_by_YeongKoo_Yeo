function gf = gfree(T,cname) 
% Estimation of Gibbs free energy of gases (kcal/gmol) at low pressure 
% input 
% T: temperature(K) (scalar or row vector) 
% cname: name or chemical formula of the compound 
% output 
% gf: Gibbs free energy of formation (kcal/gmol) 
% correlation constants (A, B) 
cv = [ 0.0 0.0; % missing component(F2)    
0.0 0.0; % missing component(Cl2)    
-71.9000 0.2500;  -26.5000 -21.3000;  -94.2000 -0.4200;  -22.3000 -1.7200;    
-12.3000 27.2000;  -58.6000 12.7000;  -33.2000 26.4000;      
0.0 0.0; % missing component(H2)    
0.0 0.0; % missing component(N2)    
0.0 0.0; % missing component(O2)    
0.0 0.0; % missing component(C2H4)     
-20.1000 24.9000;  -23.3000 49.7000;  -28.8000 74.7000;  16.7000 45.9000;    
7.8000 68.7000;  18.5000 70.6000;  -25.0000 56.5000;  10.3000 47.7000;    
-35.1000 139.0000;  24.2000 38.6000;  -50.2000 36.5000;    
-24.8000 26.9000;  -22.6000 32.1000]; 
wv = [cv(:,1) cv(:,2)*1e-3]; 
ind = compID(cname); % identification number of the component  
% computes Gibbs free energy 
if ind ~= 3, gf = wv(ind,1) + wv(ind,2)*T; 
else % ind = 3: sulfur dioxide    
gf = [];    
for j = 1:length(T)        
if T(j) >= 298 && T(j) <= 717 % temperature range: 298K~717K            
ind = 1; gfv = wv(ind,1) + wv(ind,2)*T(j);        
else % temperature range: 717K~1500K            
wv(ind,:) = [-86.8 17.7*1e-3]; gfv = wv(ind,1) + wv(ind,2)*T(j); 
end        
gf = [gf gfv];    
end 
end 
end 
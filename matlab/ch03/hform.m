function hf = hform(T,cname)
% Estimation of heat of formation of gases (kcal/gmol) at low temperature 
% input 
% T: temperature(C) (scalar or row vector) 
% cname: name or chemical formula of the compound 
% output 
% hf: heat of formation (kcal/gmol) 
% correlation constants (A, B, C) 
cv = [0.0 0.0 0.0; % missing component(F2)    
0.0 0.0 0.0; % missing component(Cl2)    
-69.6000 -5.2900 0;  -26.5000 0.8400 -1.0500;  -93.9000 -0.4300 0;  -21.9000 -0.6100 0;    
-9.3400 -6.2000 2.3600;  -57.4000 -1.7900 0;  -31.8000 -3.0400 1.1900;    
0.0 0.0 0.0; % missing component(H2)    
0.0 0.0 0.0; % missing component(N2)    
0.0 0.0 0.0; % missing component(O2)    
0.0 0.0 0.0; % missing component(C2H4)     
-15.4000 -9.5900 3.5000;  -16.4000 -14.8000 6.1300;  -20.0000 -19.1000 8.1500;    
23.7000 -15.3000 6.2700;  16.8000 -19.0000 7.8400;  25.2000 -17.6000 8.9800;    
-19.3000 -14.6000 7.1800;  16.9000 -16.7000 7.9000;  -21.6000 -32.0000 15.8000;    
29.0000 -10.1000 4.1300;  -44.9600 -11.9000 4.9800;  -24.7000 0.0334 0;    
-24.0000 2.4200 0]; 
wv = [cv(:,1) cv(:,2)*1e-3 cv(:,3)*1e-6]; 
ind = compID(cname); % identification number of the component 
% computes heat of formation 
T0 = 273.15; T = T + T0; 
if ind ~= 3    
hf = wv(ind,1) + wv(ind,2)*T + wv(ind,3)*T.^2; 
else % ind = 3: sulfur dioxide    
hf = [];    
for j = 1:length(T)        
if T(j) >= 298 && T(j) <= 717 % temperature range: 298K~717K            
ind = 1; hfv = wv(ind,1) + wv(ind,2)*T(j) + wv(ind,3)*T(j).^2;        
else % temperature range: 717K~1500K            
wv(ind,:) = [-86.9 0.32*1e-3 0]; hfv = wv(ind,1) + wv(ind,2)*T(j) + wv(ind,3)*T(j).^2;        
end        
hf = [hf hfv];    
end 
end 
end
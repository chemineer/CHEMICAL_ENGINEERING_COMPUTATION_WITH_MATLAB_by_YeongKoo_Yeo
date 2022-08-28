function sg = surtL(T,cname) 
% Calculates liquid surface tension (dyne/cm) 
% input 
% T: temperature(C) (scalar or a row vector) 
% cname: common name or chemical formula of the compound 
% output 
% sg: ratio of surface tension at T to that at T (sg = sigma/sigma1)  

% data (T1, Tc) and correlation constant (r) 
wv = [ -200.0000 -129.0000 0.8811;  20.0000 144.0000 1.0508;  30.0000 157.6000 1.1768;  
-193.0000 -140.1000 1.1441;  20.0000 31.1000 1.3015;  -93.0000 51.5000 1.0972;  
-45.0000 132.4000 1.1548;  25.0000 374.2000 0.8105;  18.2000 455.0000 0.9141;  
-256.0000 -240.2000 1.1012;  -203.0000 -146.8000 1.2123;  -202.0000 -118.5000 1.1933;  
-120.0000 9.9000 1.2760;  -168.1600 -82.6000 1.3941;  -120.0000 32.3000 1.2060;  
-90.0000 96.7000 1.1982;  20.0000 288.9400 1.2243;  20.0000 318.8000 1.2364;  
20.0000 426.0000 1.1022;  60.0000 420.0000 1.0725;  10.0000 124.9000 1.3201;  
20.0000 280.3000 1.4246;  40.0000 152.0000 1.2055;  20.0000 239.4000 0.8115;  
25.0000 263.4000 1.1824;  30.0000 283.2000 1.2278]; 
ind = compID(cname); % identification number of the component 
% compute ratio of surface tension 
T0 = 273.15; T = T + T0; T1 = wv(ind,1) + T0; Tc = wv(ind,2) + T0; 
if ind ~= 8 
sg = ((Tc - T)./(Tc-T1)).^(wv(ind,3)); 
else % ind=8: water 
sg = []; 
for j = 1:length(T)                  
if T(j) >= T0 && T(j)-T0 <= 100+T0 % 0<=T<=100 (C)                              
sg1 = 71.97;                                         
sgv = sg1*((wv(ind,2) + T0 - T(j))./(wv(ind,2) +T0-(wv(ind,1) + T0))).^(wv(ind,3));        
else % 100 <= T(C) <= 374.2            
wv(ind,:) = [100 374.2 1.169]; sg1 = 58.91;            
sgv = sg1*((wv(ind,2) + T0 - T(j))./(wv(ind,2) +T0-(wv(ind,1) + T0))).^(wv(ind,3));            
end            
sg = [sg sgv];        
end    
end 
fprintf('Ratio of surface tension(sigma/sigma1) at T = %g C is %g\n', T-T0, sg); 
end
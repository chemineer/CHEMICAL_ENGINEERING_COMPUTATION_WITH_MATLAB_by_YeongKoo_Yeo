function mu = visG(T,cname) 
% Calculates viscosity of a gas (cP) 
% input 
% T: temperature(K) (scalar or a row vector) 
% cname: common name or chemical formula of the compound 
% output 
% mu: the viscosity of the gas (cP) 
% correlation constants (A, B, C) 
cv = [22.0900 76.9000 -211.6000; 5.1750 45.6900  -88.5400;   
-3.7930 46.4500 -72.7600;  32.2800 47.4700  -96.4800;   
25.4500 45.4900 -86.4900;  -9.5540 54.4500  -96.5600;   
-9.3720 38.9900 -44.0500;  -31.8900 41.4500 -8.2720;     
5.3810 28.9800 38.4000; 21.8700 22.2000  -37.5100;      
30.4300 49.8900 -109.3000;  18.1100 66.3200 -187.9000;      
3.5860 35.1300 -80.5500;  15.9600 34.3900  -81.4000;      
5.5760 30.6400 -53.0700;  4.9120 27.1200  -38.0600;   
-15.7600 32.4500 -72.3200;  -8.4210 27.1100  -40.1800;   
-14.9800 29.0300 -1.1160;  -16.4100 32.0000 0;   
-7.7870 34.7800 -81.3000;  -4.7050 26.3200 -44.1000;   
-10.6700 34.3200 -80.8000;  -5.6360 34.4500 -3.3400;   
-6.6880 37.2600 -50.8700;  5.6980 32.7300 -40.2800]; 

wv = [cv(:,1) cv(:,2)*1e-2 cv(:,3)*1e-6]; 
ind = compID(cname); % identification number of the component 
% compute viscosity 
pn = wv(ind,1) + wv(ind,2).*T + wv(ind,3).*T.^2; mu = pn*1e-4; 
fprintf('Viscosity = %g cP\n', mu) 
end 
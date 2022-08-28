function k = condG(T,cname) 
% Calculates thermal conductivities of gases (microcal/cm/sec/K)  
% input 
% T: temperature(K) (scalar or a row vector) 
% cname: common name or chemical formula of the compound 
% output 
% k: bas thermal conductivity (microcal/cm/sec/K)  

% correlation constants (A, B, C, D) 
cv = [1.8654 19.7900 1.2400 -17.7700; 3.2500 5.8000 0.2100 -1.2500;  
-19.3100 15.1500 -0.3300 0.5500; 1.2100 21.7900 -0.8416 1.9580;  
-17.2300 19.1400 0.1308 -2.5140; -0.2600 12.6700 -0.2500 0.1600;     
0.9100 12.8700 2.9300 -8.6800; 17.5300 -2.4200 4.3000 -21.7300;  
-21.0700 16.9700 0.1700 -1.5600; 19.3400 159.7400 -9.9300 37.2900;   
0.9359 23.4400 -1.2100 3.5910; -0.7816 23.8000 -0.8939 2.3240;  
-42.0400 28.6500 0.7963 -3.2620; -4.4630 20.8400 2.8150 -8.6310;  
-75.8000 52.5700 -4.5930 39.7400; 4.4380 -1.1220 5.1980 -20.0800;  
-20.1900 8.6400 2.3400 -9.6900; 18.1400 -9.5700 5.6600 -22.2200;  
-26.3900 11.8900 1.5500 -4.3000; -31.8700 15.2600 1.7400 -4.4000;  
-20.4600 9.7400 3.7700 -16.2800; -20.5700 4.4500 4.0700 -17.3100;  
-67.9200 29.9600 1.7400 -12.2000; -18.6200 9.9500 2.9000 -12.3800;  
-5.7320 6.2900 0.5904 -3.3520; -0.4161 4.0670 0.6115 -3.5660]; 
wv = [cv(:,1) cv(:,2)*1e-2 cv(:,3)*1e-4 cv(:,4)*1e-8]; 
ind = compID(cname); % % identification number of the component 
% Computes gas thermal conductivity 
k = wv(ind,1) + wv(ind,2)*T + wv(ind,3)*T.^2 + wv(ind,4)*T.^3; 
fprintf(' Thermal conductivity = %g microcal/cm/sec/K \n', k) 
end 
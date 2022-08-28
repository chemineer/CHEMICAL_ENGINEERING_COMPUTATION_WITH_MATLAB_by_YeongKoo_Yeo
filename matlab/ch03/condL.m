function k = condL(T,cname) 
% Calculates thermal conductivities of gases (microcal/cm/sec/K)  
% input 
% T: temperature(K) (scalar or a row vector) 
% cname: common name or chemical formula of the compound 
% output 
% k: liquid thermal conductivity(microcal/cm/sec/C)  

% correlation constants (A, B, C) 
cv = 1000*[0.6215 -0.1623 -0.1184; 0.5590 -0.0483 -0.0152    
2.1408 -0.7837 0.0714;  0.4755 0.0033 -0.2143;    
0.9721 -0.2015 -0.0230;  1.0717 -0.0184 -0.0658;    
2.5513 -0.3766 -0.0294;  -0.9166 1.2547 -0.1521;    
-0.4666 0.8059 -0.0876;  -0.0201 2.4737 -5.3473;    
0.6280 -0.3689 -0.0226;  0.5838 -0.2105 -0.0483;    
0.8515 -0.2289 -0.0047;  0.7227 -0.1444 -0.0764;    
0.6993 -0.1659 -0.0049;  0.6235 -0.1268 -0.0021;    
0.4243 0.0011 -0.0090;  0.4851 -0.0538 -0.0006;    
0.5375 -0.0304 -0.0015;  0.4409 0.0867 -0.0070;    
0.3968 -0.0421 -0.0067;  0.3883 -0.0227 -0.0033;    
0.7183 -0.1872 0.0117;  0.7701 -0.1142 0.0028;    
0.3902 -0.0206 -0.0051]; 
wv = [cv(:,1) cv(:,2)*1e-2 cv(:,3)*1e-4]; 
ind = compID(cname); % identification number of the component
% Computes liquid thermal conductivity 
k = wv(ind,1) + wv(ind,2).*T + wv(ind,3).*T.^2; 
fprintf(' Thermal conductivity = %g microcal/cm/sec/C \n', k) 
end 
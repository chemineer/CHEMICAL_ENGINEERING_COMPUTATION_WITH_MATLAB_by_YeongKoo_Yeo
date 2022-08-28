function hcp = hcapL(T,cname) 
% Calculates liquid heat capacity (cal/g/C) 
% input 
% T: temperature(K) (scalar or a row vector) 
% cname: common name or chemical formula of the compound 
% output 
% hcp: liquid heat capacity(cal/g/C)    

% correlation constants (A, B, C, D) 
cv = 1.0e+04 * [-0.0000288 0.02528  -0.0331 0.1464  
-0.000013220  0.0004720  -0.0020370 0.002894000  
-0.000057370  0.0010340  -0.0040280 0.005285000  
0.000056450  0.0004798  -0.0143700 0.091195000  
-0.001930000  0.0254600  -0.1095500 0.157330000  
-0.000011210  0.0007048  -0.0035310 0.006621000  
-0.000192300  0.0031100  -0.0110900 0.013760000  
0.000067410  0.0002825  -0.0008371 0.000860100  
0.000044400  0.0001199  -0.0002738 0.000261500  
0.000379000  -0.0329800  1.2170900  -0.243480000  
-0.000106400  0.0059470  -0.0768700 0.335730000  
-0.000045870  0.0032340  -0.0395100 0.157570000  
-0.000034020  0.0006218  -0.0050120 0.012630000  
0.000123000  -0.0010330  0.0072000  -0.010730000  
0.000013880  0.0008481  -0.0056540 0.012610000  
0.000033260  0.0002332  -0.0013360 0.003016000  
-0.000148100  0.0015460  -0.0043700 0.004409000  
-0.000014610  0.0004584  -0.0013460 0.001425000  
0.000014070  0.0002467  -0.0006085 0.000592700  
-0.000068960  0.0008218  -0.0018420 0.001447000  
-0.000002618  0.0006913  -0.0034770 0.005990000  
-0.000128400  0.0013390  -0.0035100 0.003227000  
0.000037850  0.0001049  -0.0005761 0.001374000 
0.000083820  -0.0003231  0.0008296  -0.000016890  
-0.000009154  0.0003149  -0.0010640 0.001240000  
-0.000001228  0.0002058  -0.0007040 0.000861000]; 
wv = [cv(:,1) cv(:,2)*1e-3 cv(:,3)*1e-6 cv(:,4)*1e-9]; 
ind = compID(cname); % identification number of the component
% computes liquid heat capacity 
hcp = wv(ind,1) + wv(ind,2).*T + wv(ind,3).*T.^2 + wv(ind,4).*T.^3; 
fprintf(' Heat capacity = %g cal/g/C \n', hcp) 
end 
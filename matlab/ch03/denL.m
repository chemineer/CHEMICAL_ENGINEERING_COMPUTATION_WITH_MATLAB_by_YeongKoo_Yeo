function rhoL = denL(T,cname) 
% density of saturated liquid 
% input 
% T: temperature(C) (scalar or row vector) 
% Tc: critical temperature(C) 
% cname: name or chemical formula of the compound 
% output 
% rhoL: density of saturated liquid(g/cm^3) 
% coefficients A and B of the density correlation equation 
coefAB = [0.5649 0.2828 -129.0000; 0.5615 0.2720 144.0000;
0.5164 0.2554 157.6000; 0.2931 0.2706 -140.1000;        
0.4576 0.2590 31.1000; 0.4183 0.2619 51.5000;        
0.2312 0.2471 132.4000; 0.3471 0.2740 374.2000;        
0.0 0.0 0.0 % missing component(H2O2) 
0.0315 0.3473 -240.2000; 0.3026 0.2763 -146.8000;        
0.4227 0.2797 -118.5000; 0.2118 0.2784 9.9000;        
0.1611 0.2877 -82.6000; 0.2202 0.3041 32.3000;        
0.2204 0.2753 96.7000; 0.3051 0.2714 288.9400;        
0.2883 0.2624 318.8000; 0.3392 0.2761 426.0000;        
0.4094 0.3246 420.0000; 0.2614 0.2826 124.9000;        
0.2729 0.2727 280.3000; 0.2444 0.2710 152.0000;        
0.2928 0.2760 239.4000; 0.5165 0.2666 263.4000;        
0.5591 0.2736 283.2000]; 
ind = compID(cname); % identification number of the component 
% compute density 
Tr = T./coefAB(ind,3); epn = -(1 - Tr).^(2/7); rhoL = coefAB(ind,1) * coefAB(ind,2).^epn; 
fprintf('Density = %g g/cm^3\n', rhoL) 
end 
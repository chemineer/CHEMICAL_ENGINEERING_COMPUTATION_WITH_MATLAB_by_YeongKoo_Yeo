% multcstr.m 
clear all; 
ka = 10; kc = 15; V = 2500; v0 = 100; Ca0 = 2; Cb0 = 2; C0 = [Ca0 Cb0 0 0]; 
[C fval] = fsolve(@cstrmult,C0,[],ka,kc,Ca0,Cb0,v0,V); 
Ca = C(:,1); Cb = C(:,2); Cc = C(:,3); Cd = C(:,4);  
n = length(Cd); Scd = zeros(1,n); 
for i = 1:n    
if Cd(i) <= 1e-4, Scd(i) = 0; else, Scd(i) = Cc(i)/Cd(i); end 
end 
fprintf('The final concentration of each species: \n'); 
fprintf('Caf=%g, Cbf=%g, Ccf=%g, Cdf=%g\n',Ca(end),Cb(end),Cc(end),Cd(end)); 
fprintf('The final selectivity: Scdf = %g\n',Scd(end)); 
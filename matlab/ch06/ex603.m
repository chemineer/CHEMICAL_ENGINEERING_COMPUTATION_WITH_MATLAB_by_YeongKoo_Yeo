% rxnrate.m: determination of reaction order 
t = [0 2 4 7 9 18]; Ca = [1.48 1.01 0.67 0.58 0.51 0.32]; % data 
M = [t' ones(length(t),1)]; Mi = inv(M'*M)*M'; 
Y0 = Ca'; Y1 = (log(Ca))'; Y2 = (1./Ca)'; X0 = Mi*Y0; X1 = Mi*Y1; X2 = Mi*Y2; 
k0 = -X0(1); Ca0 = X0(2); k1 = -X1(1); Ca1 = exp(X1(2)); k2 = X2(1); Ca2 = 1/X2(2); 
fprintf('0th order: k = %g, Ca0 = %g\n1st order: k = %g, Ca0 = %g \n',k0,Ca0,k1,Ca1) 
fprintf('2nd order: k = %g, Ca0 = %g\n',k2,Ca2) 
% plot of each reaction rate and data 
tv = 0:0.1:20; % time interval for plotting 
Caz = -k0*tv + Ca0; Caf = exp(-k1*tv + log(Ca1)); Cas = 1./(k2*tv + 1./Ca2); 
plot(tv,Caz,':',tv,Caf,'--',tv,Cas,t,Ca,'o'), xlabel('t(sec)'), ylabel('C_A(mol/l)') 
legend('0th order','1st order','2nd order','Data','location','best') 
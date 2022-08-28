% compflr.m 
clear all; 
ka = 100; kc = 1500; Ct0 = 0.2; Ft0 = 20; alpha = 0.0019;  % data 
x0 = [10 10 0 0 1]; Wv = [0 1000]; [W x] = ode45(@pbrmult,Wv,x0, [],ka,kc,Ct0,Ft0,alpha); 
Fa = x(:,1); Fb = x(:,2); Fc = x(:,3); Fd = x(:,4); y = x(:,5); 
n = length(W); Scd = zeros(1,n); 
for i = 1:n    
if Fd(i) <= 1e-4, Scd(i) = 0; else, Scd(i) = Fc(i)/Fd(i); end 
end 
subplot(1,2,1), plot(W,Fa,W,Fb,':',W,Fc,'.-',W,Fd,'--',W,y,'.') 
xlabel('W'), ylabel('F_i'), legend('F_A','F_B','F_C','F_D','y') 
subplot(1,2,2), plot(W,Scd), xlabel('W'), ylabel('S_{C/D}') 
fprintf('Final molar flow rate of each species: \n'); 
fprintf(' Faf=%g, Fbf=%g, Fcf=%g, Fdf=%g\n',Fa(end),Fb(end),Fc(end),Fd(end)); 
fprintf('Final value of y: yf = %g\n',y(end)); fprintf('Selectivity: Scdf = %g \n',Scd(end)); 
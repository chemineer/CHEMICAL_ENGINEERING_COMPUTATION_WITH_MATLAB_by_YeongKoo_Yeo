% vpchform.m: vapor pressure of chloroform 
clear all; 
t = [0.001 10:10:150]; % temp.(deg.C) 
Pv = [65.4 101.0 159.3 242.6 357.9 513.6 718.7 983.2 1317.7 1733.9 ...     
2243.7 2859.2 3593.0 4457.9 5466.4 6631.0];  % vapor pressure (mmHg) 
n = length(t); Y = (log(Pv))'; X = [ones(n,1) (1./t)' (log(Pv)./t)']; 
K = inv(X'*X)*X'*Y; A = K(1), C = -K(3), B = A*C - K(2) 
ti = 0:0.1:150; Pvi = exp(A - B./(ti + C));  % Generate points for plot 
plot(ti,Pvi,t,Pv,'o'), xlabel('t(deg.C)'), ylabel('P_v(mmHg)') 
legend('Antoine eqn.','Data','location','best') 
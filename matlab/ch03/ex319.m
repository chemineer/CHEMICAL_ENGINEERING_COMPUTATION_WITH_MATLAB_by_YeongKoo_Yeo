% compVPeqn.m: comparison of vapor equations 
clear all; clc; 
T = [-15 -4.3 7.5 20.7 29.1 40.7 58.1 78.0 99.2] + 273.15; 
Pv = [0.667 1.333 2.666 5.333 8.000 13.33 26.66 53.33 101.32]; 
Xa = [ones(1,length(T)); 1./T; log(Pv)./T]'; ba = inv(Xa'*Xa)*Xa'*log(Pv)'; 
Aa = ba(1); Ca = -ba(3); Ba = Aa*Ca-ba(2); % Antoine eq. 
Xr = [ones(1,length(T)); 1./T; log(T); T.^6]'; br = inv(Xr'*Xr)*Xr'*log(Pv)'; 
Ar = br(1); Br = br(2); Cr = br(3); Dr = br(4); % Riedel eq. 
Xh = [ones(1,length(T)); 1./T; log(T); Pv./T.^2]'; bh = inv(Xh'*Xh)*Xh'*log(Pv)'; 
Ah = bh(1); Bh = bh(2); Ch = bh(3); Dh = bh(4); % Harlecher-Braun eq 
Ti = T(1):T(end); Pa = exp(Aa - Ba./(Ti+Ca)); 
Pr = exp(Ar + Br./Ti + Cr*log(Ti) + Dr*Ti.^6); Ph = []; 
for k = 1:length(Ti)    
P0 = 100/length(Ti) * k;    
Phf = @(x) Ah + Bh/Ti(k) + Ch*log(Ti(k)) + Dh*x/Ti(k)^2 - log(x);    
Phv = fzero(Phf,P0); Ph = [Ph Phv]; 
end 
plot(Ti,Pa,Ti,Pr,':',Ti,Ph,'.-',T,Pv,'o'), xlabel('T(C)'), ylabel('Pv(mmHg)') 
legend('Antoine','Riedel','Harlecher-Braun','Data','location','best')
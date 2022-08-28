function f = dpf(t,y,P) 
% t(1)=x1, t(2)=x2, t(3)=T 
% Required relations 
P1 = 10^(7.62231 - 1417.9/(191.15 + t(3))); P2 = 10^(8.10765 - 1750.29/(235 + t(3))); 
gam1 = 10^(1.7*t(2)^2/((2.43*t(1) + t(2))^2)); gam2 = 10^(0.7*t(1)^2/((t(1) + 0.412*t(2))^2));
k1 = gam1*P1/P; k2 = gam2*P2/P;  
% Nonlinear equations 
f = [t(1) - y(1)/k1; t(2) - y(2)/k2; t(1) + t(2) - 1];  
end 
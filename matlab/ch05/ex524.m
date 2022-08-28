% gasrelease.m: gas release from a cylinder 
T = 60; sg = 0.42; d = 0.0779; K = 2.848; % total K  
p0 = 101.3; p1 = p0 + [800 150]; n = length(p1); % pressure (kPa) 
dp1c = 1/(0.9953 + 0.9054/sqrt(K) + 0.1173/K - 0.0195/K^1.5); % limiting dp/p1 
dp1e = (p1-p0)./p1; % estimated dp/p1 
Y = 0.0415*log(K) + 0.6097; % limiting expansion factor 
for k = 1:n    
if dp1c < dp1e(k) % flow is sonic        
dp = dp1c*p1(k); % limiting pressure differnece (kPa)
Q = 3600*53.64*Y*d^2 * sqrt(dp*p1(k)/(K*(T+273.15)*sg)); % flow rate (m^3/hr)        
fprintf('Flow rate (cylinder pressure: %g kPaG) = %g m^3/sec \n',p1(k)-p0,Q);    
else % flow is subsonic        
m = (1-Y)/dp1c; % slope        
Y = 1 - m*dp1e(k); % unchoked flow        
dp = dp1c*p1(k); % limiting pressure differnece (kPa)        
Q = 3600*53.64*Y*d^2 * sqrt(dp*p1(k)/(K*(T+273.15)*sg)); % flow rate (m^3/hr)        
fprintf('Flow rate (cylinder pressure: %g kPaG) = %g m^3/sec \n',p1(k)-p0,Q);    
end 
end 
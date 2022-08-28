function dC = pfrdiff(t,C,pf) 
Pe = pf.Pe; Da = pf.Da; n = pf.n; % Retrieve data 
h = 1/n; dC = zeros(n,1); % Initialization 
% Difference equations 
for k = 1:n    
if k == 1, s = (C(k+1) - 1)/(2*h); d = (C(k+1) - 2*C(k) + 1)/(Pe*h^2);    
elseif k == n, s = 0; d = (-2*C(k) + 2*C(k-1))/(Pe*h^2);    
else, s = (C(k+1) - C(k-1))/(2*h); d = (C(k+1) - 2*C(k) + C(k-1))/(Pe*h^2); end    
dC(k) = - s + d- Da*C(k); 
end 
end 
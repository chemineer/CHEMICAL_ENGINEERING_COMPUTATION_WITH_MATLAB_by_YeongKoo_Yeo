function fD = nnDfun(D,W,rho,dP,rf,L,k,n,K) 
eD = rf./D; % roughness 
v = 4*W./(pi*rho*D.^2); % m/s 
Nre = D^n*v^(2-n)*rho/(8^(n-1)*k*((3*n+1)/(4*n))^n); % Reynolds number 
if Nre < 2100, f = 16/Nre; 
else     
f0 = 16/Nre; fe = @(x) 4*log10(Nre*x^(1-n/2))/n^0.75 - 0.4*n^(-1.2) - 1/sqrt(x);    
f = fzero(fe,f0); 
end 
Le = K*D/(4*f); fD = D - ((2*f*(L + Le))./(rho*dP) * (4*W/pi)^2).^0.2;  
end 
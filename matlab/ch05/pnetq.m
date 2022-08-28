function fun = pnetq(x,D,rho,mu,dP0,L) 
% x1=q01, x2=q12, x3=q13, x4=q23, x5=q24, x6=q34, x7=q45 
% L1=L01, L2=L12, L3=L13, L4=L23, L5=L24, L6=L34, L7=L45 
A = pi*D.^2/4; eD = 4.6e-5./D;  
for k = 1:7    
Nre = D*x(k)*rho/mu/A; % Shacham eqn.    
f = 1./(log10(eD/3.7 - (5.02./Nre).*log10(eD/3.7+14.5/Nre))).^2 /16;    
dP(k) = 32*f*rho*L(k)*(x(k))^2/(pi^2*D^5); 
end 
fun = [x(1) - x(2) - x(3);  x(2) - x(4) - x(5);  x(3) + x(4) - x(6);  x(5) + x(6) - x(7); 
dP(1) + dP(2) + dP(5) + dP(7) + dP0;  -dP(2) + dP(3) - dP(4);  dP(4) - dP(5) + dP(6)]; 
end 
% pfrgpr.m 
k = [5, 2, 10, 5]; F0 = [10, 10, 0, 0, 0, 0]; Vf = 10; Ct0 = 2; 
[V F] = ode45(@pf,[0 Vf],F0,[],k,Ct0); 
Fa = F(:,1); Fb = F(:,2); Fc = F(:,3); Fd = F(:,4); Fe = F(:,5); Ff = F(:,6); 
plot(V,Fa,'--',V,Fb,':',V,Fc,V,Fd,'.-', V,Fe,V,Ff), xlabel('V(liter)'), 
ylabel('F_i(mol/min)') 
legend('F_A','F_B','F_C','F_D', 'F_E' ,'F_F') 

function dF = pf(V,F,k,Ct0) 
% F(1) = A, F(2) = B, F(3) = C, F(4) = D, F(5) = E, F(6) = F 
% k = [k1 k2 k3 k4] 
Ft = sum(F);  
r1A = -k(1) * Ct0^3 * F(1)*F(2)^2 / (Ft^3); r2A = -k(2) * Ct0^2 * F(1)*F(2) / (Ft^2); 
r3B = -k(3)* Ct0^3 * F(2)*F(3)^2 / (Ft^3); r4C = -k(4) * Ct0^(5/3) * F(3)*F(1)^ (2/3) / (Ft^(5/3)); 
dF = [r1A + r2A + (2/3)*r4C; (5/4)*r1A + (3/4)*r2A + r3B; -r1A + 2*r3B + r4C; -(3/2)*r1A - (3/2)*r2A - r4C; -(1/2)*r2A - (5/6)*r4C; -2*r3B]; 
end 
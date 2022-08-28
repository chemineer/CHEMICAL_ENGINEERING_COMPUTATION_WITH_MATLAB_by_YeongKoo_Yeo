function fun = pipsys(x,L,d,rhg,Z,rho,mu,g,Q) 
% x(1)=Ws; x(2)=v2, x(3)=v3, x(4)=f1, x(5)=f2, x(6)=f3 
v1 = 4*Q/(pi*d(1)^2); Ws = x(1); f = [x(4) x(5) x(6)]; v = [v1 x(2) x(3)]; 
dH = f.*L.*v.^2 ./(2*g*d); Re = d.*v.*rho/mu; 
fun(1,1) = Ws + Z(1) - Z(2) - dH(1) - dH(2); 
fun(2,1) = Ws + Z(1) - Z(3) - dH(1) - dH(3); 
fun(3,1) = v(1)*d(1)^2 - v(2)*d(2)^2 - v(3)*d(3)^2; 
fun(4,1) = 1/sqrt(f(1)) + 4*log10(rhg(1)/d(1)/3.7 + 1.256/Re(1)/sqrt(f(1))); 
fun(5,1) = 1/sqrt(f(2)) + 4*log10(rhg(2)/d(2)/3.7 + 1.256/Re(2)/sqrt(f(2))); 
fun(6,1) = 1/sqrt(f(3)) + 4*log10(rhg(3)/d(3)/3.7 + 1.256/Re(3)/sqrt(f(3))); 
end 
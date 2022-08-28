function frx = pbrmult(W,x,ka,kc,Ct0,Ft0,alpha) 
% x(1)=Fa, x(2)=Fb, x(3)=Fb, x(4)=Fd, x(5)=y 
Ft = x(1) + x(2) + x(3) + x(4); 
for i = 1:4, C(i) = Ct0*x(i)*x(5)/Ft; end 
frx = [-ka*C(1)*C(2)^2 - 2*kc*C(1)^2*C(3)^3/3; -2*ka*C(1)*C(2)^2;      
ka*C(1)*C(2)^2 - kc*C(1)^2*C(3)^3; kc*C(1)^2*C(3)^3/3; - alpha*Ft/ (2*x(5)*Ft0)]; 
end 
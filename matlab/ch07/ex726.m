% btdist.m
% % z(1) = T, z(2) = L
A = [7.96681 8.04494]; B = [1668.21 1554.3]; C = [228 222.65];
c = [0.3781 0.6848]; Pt = 760; K = 5e5; x1v = [0.4 0.8]; z0 = [79 100];
[x1 z] = ode45(@bf,x1v,z0,[],Pt,A,B,C,c,K);
T = z(:,1); L = z(:,2); fprintf('Tfinal = %g, Lfinal = %g\n',T(end),L(end))
subplot(1,2,1), plot(x1,T), xlabel('x_1'), ylabel('T(C)')
subplot(1,2,2), plot(x1,L), xlabel('x_1'), ylabel('L(kgmole)')

function dz = bf(x1,z,Pt,A,B,C,c,K) 
% z(1) = T, z(2) = L 
x2 = 1-x1; x = [x1 x2]; Pj = 10.^(A - B./(z(1) + C)); 
gam(1) = 10^((1-x1)^2*(c(1) + 2*x1*(c(2)-c(1)))); gam(2) = 10^((1-x2)^2*(c(2) + 2*x2*(c(1)-c(2))));  
k = gam.*Pj/Pt; 
dz(1,1) = K*(1 - k(1)*x1 - k(2)*x2); dz(2,1) = z(2)/(x1*(k(1)-1)); 
end 

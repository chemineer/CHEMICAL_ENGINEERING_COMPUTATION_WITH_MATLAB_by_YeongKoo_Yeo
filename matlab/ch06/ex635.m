% mbrmult.m
clear all;
k1 = 2; k2 = 3; Ct0 = 0.8; Fb0 = 4; Vt = 50; Fai = 4; % data
Rb = Fb0/Vt; x0 = [Fai 0 0 0]; Vspan = [0 Vt];
[V x] = ode45(@mbf,Vspan,x0,[],k1,k2,Ct0,Rb);
Fa = x(:,1); Fb = x(:,2); Fd = x(:,3); Fu = x(:,4); n = length(V); Sdu = zeros(1,n);
for i = 1:n
if Fu(i) <= 1e-6, Sdu(i) = 0; else, Sdu(i) = Fd(i)/Fu(i); end
end
subplot(1,2,1), plot(V,Fa,V,Fb,':',V,Fd,'.-',V,Fu,'--'), xlabel('V(dm^3)'), ylabel('F_i(mol/s)')
legend('F_A','F_B','F_D','F_U','Location','best')
subplot(1,2,2), plot(V,Sdu), xlabel('V(dm^3)'), ylabel('S_{D/U}')
fprintf('Final molar flow rates: \n');
fprintf('Faf = %g,  Fbf = %g,  Fdf = %g,  Fuf = %g\n',Fa(end),Fb(end),Fd(end),Fu(end));
fprintf('Overall selectivity: Sduf = %g\n',Sdu(end));
function frx = mbf (V,x,k1,k2,Ct0,Rb)
% x(1)=Fa, x(2)=Fb, x(3)=Fd, x(4)=Fu
C = Ct0*x/sum(x);
frx = [-k1*C(1)^2*C(2) - k2*C(1)*C(2)^2; -k1*C(1)^2*C(2) - k2*C(1)*C(2)^2 + Rb;
k1*C(1)^2*C(2); k2*C(1)*C(2)^2];
end

% bpbin.m: bubble point P for binary system 
T = 60+273.15; R = 8.314; % t=60 deg.C 
P1s = 83.25e3; P2s = 37.97e3; A12 = 0.59; A21 = 1.42;  
B11 = -963e-6; B22 = -1523e-6; B12 = 52e-6; % (m^3/mol) 
x = [0:0.01:1]; n = length(x); % x1 
gam1 = exp((1-x).^2 .*(A12 + 2*(A21 - A12)*x)); % activity coefficient gamma1 
gam2 = exp(x.^2 .*(A21 + 2*(A12 - A21)*(1-x))); % activity coefficient gamma2 
d12 = 2*B12 - B11 - B22; 
for k = 1:n % z(1)=y1, z(2)=P    
f = @(z) [x(k)*gam1(k)*P1s-z(1)*z(2)*exp((B11*(z(2)-P1s)+z(2)*d12*(1-z(1))^2)/R/T);...        
(1-x(k))*gam2(k)*P2s-(1-z(1))*z(2)*exp((B22*(z(2)-P2s)+z(2)*d12*(z(1))^2)/R/T)];    
z0 = [0.5, 5e4]; z = fsolve(f,z0); y(k) = z(1); P(k) = z(2); 
end 
plot(x,P,':',y,P), xlabel('x_1, y_1'), ylabel('P(Pa)') 
legend('Bubble point','Dew point','location','best') 
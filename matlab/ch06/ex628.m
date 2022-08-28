% mcatb.m 
% y1=Xa1, y2=Xa2, y3=Xa3, y4=a2, y5=a3 
clear all; 
u=8; Ap=12; k=30; Ka=5; Ca0=0.2; kds=17.5; kdp=140; 
zspan = [0 6]; y0 = [0 0 0 1 1]; 
dydz = @(z,y) [(1/(1+Ap*sqrt(z/u)))*k*Ca0*(1-y(1))/(1+Ka*Ca0*(1-y(1)))/u;                
y(4)*k*Ca0*(1-y(2))/(1+Ka*Ca0*(1-y(2)))/u;                
y(5)*k*Ca0*(1-y(3))/(1+Ka*Ca0*(1-y(3)))/u;                
-kds*y(4)^2/u;                 
-kdp*Ca0*y(3)*y(5)/u]; 
[z,y] = ode45(dydz, zspan, y0); 
subplot(1,2,1), plot(z,y(:,1),'-',z,y(:,2),':',z,y(:,3),'--');  
xlabel('z(m)'), ylabel('X_A'), legend('X_{A1}','X_{A2}','X_{A3}'); 
subplot(1,2,2), plot(z,y(:,4),'-',z,y(:,5),':'), xlabel('z(m)'), ylabel('a'), legend('a_s','a_p'); 
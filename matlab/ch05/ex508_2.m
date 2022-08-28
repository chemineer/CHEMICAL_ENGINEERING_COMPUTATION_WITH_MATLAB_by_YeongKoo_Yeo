% vwfig.m: calculates flow velocities and flow rates 
clear all; T = 60; rf = 0.00015; dz = 300; dP = -150; 
D = [4.026 5.047 6.065 7.981]; nD = length(D);  
L = 500:500:10000; nL = length(L); vr = []; v0 = 10; qr = []; 
for i = 1:nD    
for j = 1:nL         
v(j) = fzero(@vwfun, v0, [], T,L(j),D(i),rf,dz,dP);         
q(j) = (7.481*60)*(pi*v(j).*(D(i)/12).^2)/4; % gpm    
end    
vr = [vr v']; qr = [qr q'];  
end 
figure(1), plot(L,vr(:,1),'o',L,vr(:,2),'*',L,vr(:,3),'+',L,vr(:,4),'d') 
legend('D=4in','D=5in','D=6in','D=8in'), axis([500 10000 0 20]), xlabel('L (ft)'), ylabel('Velocity, v(ft/s)') 
figure(2), plot(L,qr(:,1),'o',L,qr(:,2),'*',L,qr(:,3),'+',L,qr(:,4),'d') 
legend('D=4in','D=5in','D=6in','D=8in'), xlabel('L(ft)'), ylabel('Flow rate, q(gpm)') 
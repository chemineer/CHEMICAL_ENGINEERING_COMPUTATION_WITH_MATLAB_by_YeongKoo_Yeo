% unsflow.m: unsteady flow through pipe 
clear all; 
x = linspace(0,1,100); % 100 values from 0 to 1 
t = linspace(0,1,60); % time span: 60 sec 
m = 1; res = pdepe(m,@fleqn,@itcon,@bncon,x,t); u = res(:,:,1); 
subplot(1,2,1), surf(x,t,u), colormap(gray), xlabel('\xi'), ylabel('\tau'), 
zlabel('\phi') 
shading interp, hold on  % surface plot 
surf(-x,t,u), colormap(gray), xlabel('\xi'), ylabel('\tau'), zlabel('\phi') 
shading interp, hold off  % axisymmetric plot 
subplot(1,2,2) 
% plot tau as a function of xi 
for k = 1:length(t), plot(x,u(k,:),'k'), xlabel('\xi'), ylabel('\tau'), hold 
on; end 
for k = 1:length(t), plot(-x,u(k,:),'k'), xlabel('\xi'), ylabel('\tau'), hold 
on; end 
hold off 
function [c,f,s] = fleqn(x,t,u,DuDx) 
c = (1); f = DuDx; s = 4; 
end 
function u0 = itcon(x) 
u0 = 0; 
end 
function [pl,ql,pr,qr] = bncon(xl,ul,xr,ur,t) 
pl = 0; ql = 1; % 1st boundary 
pr = ur; qr = 0; % 2nd boundary 
end 
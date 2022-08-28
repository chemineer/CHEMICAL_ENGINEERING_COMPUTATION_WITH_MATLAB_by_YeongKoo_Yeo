% bzacP.m: vapor phase composition and activity coefficient 
% for benzene/acetic acid system 
x1 = [0.0 0.0069 0.1565 0.3396 0.4666 0.6004 0.7021 0.8286 0.8862 0.9165 0.9561 0.9840 1.0]; 
P = [57.52 58.2 126.00 175.30 189.50 224.30 236.00 250.20 259.00 261.11 264.45 266.53 271.00]; 
c = polyfit(x1,P,4) % regression polynomial (4th order) 
dc = polyder(c) % differentiation of the polynomial 
x10 = 1e-5; xinv = [x10 1]; y10 = x10*dc(end)/c(end); 
[x,y] = ode45(@bzacfun,xinv,y10); % solve the differential eqn. 
P1v = P(end); P2v = P(1); Px = polyval(c,x); 
gam1 = y.*Px./x/P1v; gam2 = (1-y).*Px./(1-x)/P2v; 
subplot(1,2,1), plot(x,y), xlabel('x_1'), ylabel('y_1') % x1-y1 graph 
subplot(1,2,2), plot(x,gam1,x,gam2,'.-'), xlabel('x_1'), ylabel('\gamma') 
legend('\gamma_1','\gamma_2','Location','best'), axis([0 1 0 4])
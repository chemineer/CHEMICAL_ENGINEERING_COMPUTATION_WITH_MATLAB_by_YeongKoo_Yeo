% catwtx.m 
clear all; 
k = 0.00392; Fa0 = 0.1362; Fb0 = 0.068; P0 = 10; alpha = 0.0367; % data 
Fi = Fb0*(79/21); Ft0 = Fa0 + Fb0 + Fi; ya0 = Fa0/Ft0; 
Pa0 = ya0*P0; kp = k*Pa0*(0.5)^(2/3); delta = -0.5; epsilon = ya0*delta; 
Wfa = 15; Wfb = 25; % initial guess of catalyst weight 
z0 = [0 1]; Xf = 0.6; 
while (abs(Wfa-Wfb) >= 1e-3)     
Wfm = (Wfa+Wfb)/2; [Wa Za] = ode45(@pbconv,[0 Wfa],z0, [],kp,Fa0,epsilon,alpha); 
[Wm Zm] = ode45(@pbconv,[0 Wfm],z0,[],kp,Fa0,epsilon,alpha); 
[Wb Zb] = ode45(@pbconv,[0 Wfb],z0,[],kp,Fa0,epsilon,alpha); 
Xa = Za(end,1); Xm = Zm(end,1); Xb = Zb(end,1); 
if (Xa-Xf)*(Xm-Xf) < 0, Wfb = Wfm;  
else, Wfa = Wfm; end 
end 
X = Zm(:,1); y = Zm(:,2); W = Wm; 
fprintf('Catalyst weight = %g, Conversion = %g\n',Wfm, X(end)); 
fm = (1 + epsilon*X)./y; % ratio of the volumetric flow rate 
rp = -kp*(1 - X).*y./(1 + epsilon*X); % reaction rate 
subplot(1,2,1), plot(Wm,X,Wm,y,':',Wm,fm,'--'), legend('X','y','f'), xlabel('W(kg)'), grid on 
subplot(1,2,2), plot(Wm,-rp), xlabel('W(kg)'), ylabel('-r_A'), grid on 

function dz = pbconv(w,z,kp,Fa0,epsilon,alpha) 
% x = z(1), y = z(2) 
dz = [kp*(1-z(1))*z(2)/(Fa0*(1 + epsilon*z(1))); - alpha*(1 + epsilon*z(1))/ (2*z(2))]; 
end 
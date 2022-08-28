% encdtank.m 
Cv1 = 1.2e-2; Cv2 = Cv1; P0 = 1.014e5; P10 = 1.38e5; P30 = 1.08e5;  
P1 = P10; P3 = P30; z0 = 3.05; rho = 1e3; A = 0.465; g = 9.8; P20 = P0 + rho*g*z0; 
cv = 0.2; R = 1.987; Mw = 29; V0 = 2.83; Vg0 = 1.415;  
k1 = sqrt(P10/(P10-P20)); k2 = k1*Cv2/Cv1;  
k3 = rho*g*z0/P10; k4 = V0/Vg0; k5 = 0.735; k6 = R/(cv*Mw); 
tspan = [0 0.4]; z0 = 1; [t z] = ode45(@cltank,tspan,z0,[],k1,k2,k3,k4,k5,k6,P10,P1,P3); 
Vg = k4 - z; Tg = (1./Vg).^k6; Pg = k5*Tg./Vg; P1s = P1/P10; P3s = P3/P10; P2s = Pg + k3*z;  
F1s = k1*sqrt(P1s - P2s); F2s = k2*sqrt(P2s - P3s); 
plot(t,z,t,F1s,'--',t,F2s,'.-',t,P2s,':'), legend('z','F_1','F_2','P_2'), xlabel('t(dimensionless)') 

function dz = cltank(t,z,k1,k2,k3,k4,k5,k6,P10,P1,P3) 
P1s = P1/P10; P3s = P3/P10; Vg = k4 - z; Tg = (1/Vg)^k6; Pg = k5*Tg/Vg; 
P2s = Pg + k3*z; F1s = k1*sqrt(P1s - P2s); F2s = k2*sqrt(P2s - P3s); dz = F1s - F2s; 
end 
% adbpfr.m 
clear all;  
Ca0 = 9.3; Fa0 = 146.7; T1 = 360; T2 = 333; k1 = 31.1; K2 = 3.03; 
E = 65700; R = 8.314; dH = -6900; Vspan = [0 5]; X0 = 0; 
[V X] = ode45(@fad,Vspan,X0,[],Ca0,Fa0,T1,T2,k1,K2,E,R,dH); 
T = 330 + 43.4265*X; k = k1*exp(E*(1/T1 - 1./T)/R); Kc = K2*exp(dH*(1/T2 - 1./ T)/R); 
ra = -k*Ca0.*(1 - (1 + 1./Kc).*X); Xe = Kc./(1+Kc); 
subplot(2,2,1), plot(V,T), xlabel('V'), ylabel('T(K)') 
subplot(2,2,2), plot(V,-ra), xlabel('V'), ylabel('-r_A') 
subplot(2,2,3), plot(V,X,V,Xe,'--'), xlabel('V'), ylabel('X,X_e'), legend('X','Xe') 
fprintf('Conversion (X) and equilibrium conversion (Xe): Xf = %g, Xef = %g\n',X(end),Xe(end)); 
fprintf('Final temperature: Tf = %g\n',T(end)); 
fprintf('Final reaction rate: raf = %g\n',-ra(end)); 

function dX = fad (V,X,Ca0,Fa0,T1,T2,k1,K2,E,R,dH) 
T = 330 + 43.4265*X; k = k1*exp(E*(1/T1 - 1/T)/R); Kc = K2*exp(dH*(1/T2 - 1/ T)/R); 
ra = -k*Ca0*(1 - (1 + 1/Kc)*X); dX = -ra/Fa0; 
end 
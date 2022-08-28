% isomnbt.m 
clear all;  
Ca0 = 9.3; Fa0 = 14.67; T1 = 360; T2 = 333; k1 = 31.1; K2 = 3.03; % data 
E = 65700; R = 8.314; dH = -6900; Cp0 = 158.889; Cpc = 28; Ua = 5000; m = 500; hx = 'cn'; % data and conditions 
criT = 1e-3; Taf = 310; errT = 10; Vi = [0 5]; 
if hx == 'cn' % guess Ta(z1) 
  Ta01 = 305; Ta02 = 320; Ta0m = (Ta01+Ta02)/2;     
while errT >= criT        
z01 = [Ta01 0 305]; z02 = [Ta02 0 305]; z0m = [Ta0m 0 305];         
[V z1] = ode45(@exbco,Vi,z01, [],Ca0,Fa0,Cp0,Cpc,Ua,m,T1,T2,k1,K2,E,R,dH,hx);         
[V zm] = ode45(@exbco,Vi,z0m, [],Ca0,Fa0,Cp0,Cpc,Ua,m,T1,T2,k1,K2,E,R,dH,hx);         
[V z2] = ode45(@exbco,Vi,z02, [],Ca0,Fa0,Cp0,Cpc,Ua,m,T1,T2,k1,K2,E,R,dH,hx);        
if (z1(end,1)-Taf)*(zm(end,1)-Taf) < 0, Ta02 = Ta0m; Ta0m = (Ta01+Ta02)/2;        
else, Ta01 = Ta0m; Ta0m = (Ta01+Ta02)/2; end        
errT = abs(zm(end,1) - Taf);    
end 
else    
z0 = [310 0 305]; [V zm] = ode45(@exbco,Vi,z0, [],Ca0,Fa0,Cp0,Cpc,Ua,m,T1,T2,k1,K2,E,R,dH,hx); 
end 
Ta = zm(:,1); X = zm(:,2); T = zm(:,3); k = k1*exp(E*(1/T1 - 1./T)/R); Kc = K2*exp(dH*(1/T2 - 1./T)/R); 
ra = -k*Ca0.*(1 - (1 + 1./Kc).*X); Xe = Kc./(1+Kc); 
subplot(2,2,1), plot(V,T,V,Ta,'--'), xlabel('V'), ylabel('T(K)'), legend('T','T_a') 
subplot(2,2,2), plot(V,X,V,Xe,'--'), xlabel('V'), ylabel('X,X_e'), legend('X','Xe') 
subplot(2,2,3), plot(V,-ra), xlabel('V'), ylabel('-r_A') 
fprintf('Conversion (X) and equilibrium conversion (Xe): Xf = %g, Xef = %g\n',X(end),Xe(end)); 
fprintf('Final T and Ta: Tf = %g, Taf = %g\n',T(end), Ta(end)); 
fprintf('Final reaction rate: raf = %g\n',-ra(end)); 
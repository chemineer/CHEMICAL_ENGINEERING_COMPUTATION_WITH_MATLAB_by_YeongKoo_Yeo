% cstrun.m 
% data for CSTR 
Caf = 10; D = 2.335; Fi = 10; Tf = 25; Tj = 25; dH = 5960; A1 = 3.49308e7;  
E = 11843; rCp = 500; Ui = 70; R = 1.987; 
% initial guess and time span 
X0 = [1 5 20]; tspan = [0 20]; 
[t X] = ode45(@cstrxn,tspan,X0,[],Caf,D,Fi,Tf,Tj,dH,A1,E,rCp,Ui,R); 
% plot results 
subplot(2,2,1), plot(t,X(:,1)), xlabel('t(hr)'), ylabel('h(m)'), grid  
subplot(2,2,2), plot(t,X(:,2)), xlabel('t(hr)'), ylabel('C_A(kmol/ m^3)'), grid  
subplot(2,2,3), plot(t,X(:,3)), xlabel('t(hr)'), ylabel('T(deg.C)'), grid  

function dX = cstrxn(t,X,Caf,D,Fi,Tf,Tj,dH,A1,E,rCp,Ui,R) 
% X(1) = h, X(2) = Ca, X(3) = T(deg.C) 
Sa = (pi/4)*D^2; Sh = Sa + pi*D*X(1); 
a1 = Fi/Sa/X(1); k1 = A1*exp(-E/R/(X(3)+273.15)); 
b1 = dH/rCp; c1 = Ui*Sh/(rCp*Sa*X(1)); 
dX = [Fi/Sa - sqrt(10*X(1)/Sa); a1*(Caf - X(2)) - k1*X(2);    
a1*(Tf - X(3)) + b1*k1*X(2) - c1*(X(3) - Tj)]; 
end
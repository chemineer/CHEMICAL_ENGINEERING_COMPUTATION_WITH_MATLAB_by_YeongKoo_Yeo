% batrun.m 
X0 = [1.5 0 50]; tspan = [0 6000]; 
[t X] = ode45(@batrxn,tspan,X0); 
subplot(1,2,1), plot(t,X(:,1),t,X(:,2),'--') 
xlabel('t(s)'), ylabel('C(kmol/m^3)'), legend('C_A(t)','C_B(t)')  
subplot(1,2,2), plot(t,X(:,3)), xlabel('t(s)'), ylabel('T(deg.C)') 

function dX = batrxn(t,X) 
% X(1) = Ca, X(2) = Cb, X(3) = T(deg.C) 
batdat; % supply data 
a1 = dH1/(rho*Cp); a2 = dH2/(rho*Cp); a3 = Uj*AjV/(rho*Cp); a4 = Uc*AcV/(rho*Cp); 
k1 = A1*exp(-E1/R/(273.15+X(3))); k2 = A2*exp(-E2/R/(273.15+X(3))); 
dX = [-k1*X(1)^2; k1*X(1)^2 - k2*X(2); a1*k1*X(1)^2 + a2*k2*X(2) + a3*(Ts-X(3)) - a4*(X(3)-Tc)]; 
end 
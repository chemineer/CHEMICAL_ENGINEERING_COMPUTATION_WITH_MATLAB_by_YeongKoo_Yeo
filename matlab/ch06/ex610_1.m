% sscstr.m: steady-state solution for a CSTR 
clear all; 
Caf = 10; D = 2.335; Fi = 10; Tf = 25; Tj = 25; dH = 5960; A1 = 3.49308e7;  
E = 11843; rCp = 500; Ui = 70; R = 1.987; 
Sa = (pi/4)*D^2; hs = Fi^2 / (10*Sa); Sh = Sa + pi*D*hs;  
a1 = Fi/Sa/hs; b1 = dH/rCp; c1 = Ui*Sh/(rCp*Sa*hs); 
fs = @(X) [a1*(Caf-X(1))-A1*exp(-E/R/(X(2)+273.15))*X(1); % X(1) = Cas, X(2) = Ts    
a1*(Tf - X(2)) + b1*A1*exp(-E/R/(X(2)+273.15))*X(1) - c1*(X(2) - Tj)]; 
% Solve nonlinear equations for three different initial guesses. 
x01 = [8 20]; x1 = fsolve(fs,x01); x02 = [5 70]; x2 = fsolve(fs,x02); x03 = [2 120]; x3 = fsolve(fs,x03); 
fprintf('Initial guess: Cas = %g, Ts = %g',x01(1), x01(2)); 
fprintf('  Steady-state values: Cas = %g, Ts = %g\n',x1(1), x1(2)); 
fprintf('Initial guess: Cas = %g, Ts = %g',x02(1), x02(2)); 
fprintf('  Steady-state values: Cas = %g, Ts = %g\n',x2(1), x2(2)); 
fprintf('Initial guess: Cas = %g, Ts = %g',x03(1), x03(2)); 
fprintf('  Steady-state values: Cas = %g, Ts = %g\n',x3(1), x3(2)); 
% Generate heat and temperature profiles. 
Ts = [20:0.1:120]; k1s = A1*exp(-E/R./(Ts+273.15)); Cas = a1*Caf./(k1s+a1); 
Qr = Ui*Sh*(Ts-Tj) + Fi*rCp*(Ts-Tf); Qg = dH*Sa*hs*Cas.*k1s; 
Tjs = Ts + (rCp*Fi*(Ts-Tf) - dH*Sa*hs*Cas.*k1s)/(Ui*Sh); 
subplot(1,2,1) % Heat profile 
plot(Ts,Qr,':',Ts,Qg), xlabel('Reactor temp.(deg.C)'), ylabel('Q(kcal/h)') 
legend('Qr','Qg','location','best'), grid 
subplot(1,2,2) % Temperature profile 
plot(Tjs,Ts), axis([0 60 25 110]), grid, xlabel('T_j(deg.C)'), ylabel('T(deg.C)') 
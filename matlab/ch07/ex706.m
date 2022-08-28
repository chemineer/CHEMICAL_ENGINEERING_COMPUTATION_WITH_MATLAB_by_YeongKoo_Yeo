% crxf.m
Cas = 3e-5; Ct=4e-5; L=0.2; De=0.01; k=8e4; Kc=6e5; zspan = [0 L]; critN = 1e-10;
errN = 1;
Na1 = 2e-6; Na2 = 5e-6;
while errN > critN
Nam = (Na1+Na2)/2;
[z x1] = ode45(@crf,zspan,[Na1,Cas],[],Ct,k,Kc,De);
[z x2] = ode45(@crf,zspan,[Na2,Cas],[],Ct,k,Kc,De);
[z xm] = ode45(@crf,zspan,[Nam,Cas],[],Ct,k,Kc,De);
if x1(end,1)*xm(end,1) < 0, Na2 = Nam;
else, Na1 = Nam; end
errN = abs(Na1 - Na2);

end
Na = xm(:,1); Ca = xm(:,2); ras = Cas^2 - (Ct-Cas)/Kc; effc = Na(1)/(L*k*ras);
fprintf('Effectiveness factor: %7.5f\n', effc);
fprintf('Molar flux of A at z=0: %12.8f\n', Na(1));
subplot(1,2,1), plot(z,Ca), xlabel('z(cm)'), ylabel('C_A(gmol/cm^3)'), grid
subplot(1,2,2), plot(z,Na), xlabel('z(cm)'), ylabel('N_A(gmol/s/cm^2)'), grid

function dx = crf(z,x,Ct,k,Kc,De)
% x1=Na, x2=Ca
dx = [-k*(x(2)^2 - (Ct-x(2))/Kc); x(1)*(x(2)/Ct - 2)/(2*De)];
end

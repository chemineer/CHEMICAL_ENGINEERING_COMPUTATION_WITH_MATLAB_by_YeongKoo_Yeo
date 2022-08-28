% compintp.m: compare performance of various interpolation methods 
t  =  [0,2,5,6,13,20,24,30,35,41,50]; 
c  =  [0.86,0.61,0.47,0.39,0.25,0.18,0.15,0.12,0.10,0.09,0.08]; 
ti  =  [8,15,25,32]; 
cL  =  lagrintp(t,c,ti); % Lagrange method 
cN  =  newtintp(t,c,ti); % Newton method 
cC  =  cspintp(t,c,ti); % cubic spline method 
fprintf('Lagrange method: C  =  '); fprintf('%g ', cL); 
fprintf('\nNewton method: C  =  '); fprintf('%g ', cN); 
fprintf('\nCubic spline method: C  =  '); fprintf('%g ', cC); 
plot(t,c,'o-',ti,cL,'*',ti,cN,'<',ti,cC,'kd') 
xlabel('t(min)'), ylabel('C_A(mol/l)') 
legend('Data','Lagrange','Newton','Cubic spline','location','best')
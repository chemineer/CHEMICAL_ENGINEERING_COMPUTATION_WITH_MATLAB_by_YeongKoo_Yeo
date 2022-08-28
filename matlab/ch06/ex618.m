% pfrmxplot.m 
% Plot concentrations 
k = [55.2 30.2]; C0 = [0.021 0.0105 0]; tf = 0.5; % data and initial conditions 
[t C] = ode45(@pfrmx,[0 tf],C0,[],k,C0); Ch = C(:,1); Cm = C(:,2); Cx = C(:,3); 
plot(t,Ch,'--',t,Cm,':',t,Cx), xlabel('\tau(hr)'), ylabel('C(lbmol/ ft^3)'), legend('C_H','C_M','C_X') 
% Find optimum t 
Cxm = max(Cx); ti = find(Cx == Cxm); opmt = t(ti);  
fprintf('Optimum residence time = %g\n', opmt); 
fprintf('Maximum concentration of m-xylene = %g\n', Cxm); 

function dC = pfrmx(t,C,k,C0) 
% C(1) = Ch, C(2) = Cm, C(3) = Cx 
% k = [k1 k2], C0 = [Ch(0) Cm(0) Cx(0)] 
dC = [-k(1)*C(1)^0.5*C(2) - k(2)*C(3)*C(1)^0.5; -k(1)*C(1)^0.5*C(2);    
k(1)*C(1)^0.5*C(2) - k(2)*C(3)*C(1)^0.5]; 
end
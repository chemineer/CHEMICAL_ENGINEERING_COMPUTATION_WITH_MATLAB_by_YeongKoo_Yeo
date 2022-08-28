% compflow.m: compressible fluid flow 
d = 0.614; Lst = 20; f = 0.026; K = 2.026; Z = 0.9; T = 100; Mw = 19.5; 
P1 = 1124.7; P2 = 414.7; k = 1.27; mu = 0.012; T = T+460; % T(R) 
D = d/12; Area = pi*D^2/4; % cross sectional area (ft^2) 
rho = P1*Mw/(10.72*Z*T); % density (lb/ft^3) 
delP = P1 - P2; Kp = f*Lst/D; Kt = K + Kp; %Ktotal 
G = 1335.6*d^2 * sqrt(rho*(P1^2 - P2^2)/(Kt + 2*log(P1/P2))/P1); 
Nre = 6.31*G/(D*mu); Vg = 0.0509*G/(rho*d^2); 
Vs = 223*sqrt(k*T/Mw); Mach = Vg/Vs; Sg = Mw/29; R = 1544/(29*Sg);  
Pc = (G/(11400*d^2))*sqrt(R*T/(k*(k+1))); 
G, Vg, Mach, Pc 
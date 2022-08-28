% vppreos.m: Vp by PR EOS 
R = 83.14; w = 0.224; Tc = 304.2; Pc = 73.83; % data for CO2 
T = 288.15; Tr = T/Tc; 
a = 0.45724*R^2*Tc^2*(1+(0.37464+1.54226*w-0.26992*w^2)*(1-sqrt(Tr)))^2/Pc;  
b = 0.0778*R*Tc/Pc; 
V = 40:0.5:400; P = R*T./(V - b) - a./(V.^2 + 2*b*V - b^2); % P by Peng-Robinson EOS 
plot(V,P), axis([40 400 0 100]), grid, xlabel('V(cm^3/mol'), ylabel('P(bar)') 
% Determine Pv using bisection method 
x0l = 60; x0g = 300;  % initial guess for Vl and Vg 
xa = 43; fa = @(x) R*T./(x - b) - a./(x.^2 + 2*b*x - b^2) - xa; 
xb = 55; fb = @(x) R*T./(x - b) - a./(x.^2 + 2*b*x - b^2) - xb; 
xm = (xa + xb)/2; fm = @(x) R*T./(x - b) - a./(x.^2 + 2*b*x - b^2) - xm; 
Ppr = @(x) R*T./(x - b) - a./(x.^2 + 2*b*x - b^2); 
Vla = fsolve(fa,x0l); Vga = fsolve(fa,x0g); da = quad(Ppr,Vla,Vga) - xa*(Vga - Vla); 
Vlb = fsolve(fb,x0l); Vgb = fsolve(fb,x0g); db = quad(Ppr,Vlb,Vgb) - xb*(Vgb - Vlb); 
Vlm = fsolve(fm,x0l); Vgm = fsolve(fm,x0g); dm = quad(Ppr,Vlm,Vgm) - xm*(Vgm - Vlm); 
crit = abs(xa - xb); 
while crit > 1e-3  % Pv by bisection method    
if da*dm < 0, xb = xm; db = dm;    
elseif dm*db < 0, xa = xm; da = dm; end    
xm = (xa + xb)/2; Vlm = fsolve(fm,x0l); Vgm = fsolve(fm,x0g);     
dm = quad(Ppr,Vlm,Vgm) - xm*(Vgm - Vlm); crit = abs(xa - xb); 
end 
Pv = xm 

% The vapor pressure by the extended Antoine equation
T = 288.15; Pv = 10^(47.544-1792.2/T -16.559*log10(T) + 0.013833*T)/ 750.0615
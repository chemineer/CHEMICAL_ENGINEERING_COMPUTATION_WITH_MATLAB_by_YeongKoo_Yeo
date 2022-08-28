% twophdP.m: pressure drop in two-phase flow 
clear all;  
rhol = 66.7; rhog = 2.98; mul = 1; mug = 0.02; sigl = 70; % data 
L = 26400; d = 6.065; Wl = 77956; Wg = 12434; D = d/12; eD = 0.00015/D; 
% liquid/vapor phase friction factor 
Nrel = 6.31*Wl/(d*mul); Nreg = 6.31*Wg/(d*mug); 
if Nrel <= 2100, fL = 64/Nrel;  
else, Av = eD/3.7 + (6.7/Nrel)^0.9; fL = 4./(-4*log10(eD/3.7 - 5.02*log10(Av)/ Nrel)).^2; end 
if Nreg <= 2100, fG = 64/Nreg; 
else, Av = eD/3.7 + (6.7/Nreg)^0.9; fG = 4./(-4*log10(eD/3.7 - 5.02*log10(Av)/ Nreg)).^2; end 
% Determine two-phase flow modulus (Yg) and flow regimes 
fr = {'stratified','wave','plug','slug','bubble','annular','dispersed'}; 
rind = 7;  
delPl = 3.66e-4 * fL * Wl^2/(d^5 * rhol); delPg = 3.66e-4 * fG * Wg^2/(d^5 * rhog); 
[Yg rind] = twophmod(rhol,rhog,mul,sigl,Wl,Wg,delPl,delPg,d); 
if rind == 2            
fH = exp(0.2111*log(Wl*mul/(Wg*mug)) - 3.993);    
delPt = 3.66e-4 * fH * Wg^2 * L/(d^5 * rhog * 100); 
else, delPt = Yg*delPg*L/100; end 
fprintf('Flow regime: %s\n', fr{rind}); fprintf('Two-phase flow modulus (Yg): %g\n', Yg); 
fprintf('Total pressure drop(psi): %g\n', delPt); 

function [Yg rind] = twophmod(rhol,rhog,mul,sigl,Wl,Wg,delPl,delPg,d) 
% Calculation of two-phase flow modulus (Yg) 
% input: 
% rhol,rhog: density of liquid and vapor (lb/ft^3), mul: liquid viscosity (cP) 
% sigl: liquid surface tension (dyne/cm), Wl,Wg: flow rates of liquid and vapor (lb/h) 
% delPl,delPg: pressure drop for liquid-only/vapor-only flow, d: pipe inside diameter (in) 
% output: 
% Yg: two-phase flow modulus,  rind: flow regime index 
% Calculate Baker parameter Bx and By 
D = d/12; A = pi*D^2/4; Bx = 531*(Wl/Wg)*(sqrt(rhol*rhog)/(rhol^(2/3)))*(mul^ (1/3)/sigl); 
By = 2.16*(Wg/A)/sqrt(rhol*rhog); C = twophreg(Bx); 
% Classify flow regimes and calculate Yg for each flow regime 
x = sqrt(delPl/delPg); Wa = Wl/A; 
if By <= C(1)    
if By <= C(2), [Yg rind] = strat(x,Wa); % By < C1,C2    
else, [Yg rind] = wave(x,Wa); end % C2 < By < C1    
else % C1 < By        
if By < C(5) % C1 < By < C5            
if By < C(6), [Yg rind] = plug(x,Wa);  % C1 < By < C5,C6            
else % C1,C6 < By < C5                
if By < C(4), [Yg rind] = slug(x,Wa); % C1,C6 < By < C4,C5                
else % C1,C4,C6 < By < C5                    
if By <= C(3), [Yg rind] = annul(x,d); % C1,C4,C6 < By < C3,C5                    
else, [Yg rind] = dispr(x); % C1,C3,C4,C6 < By < C5                    
end                
end            
end        
else % C1,C5 < By            
if Bx > 150, [Yg rind] = bubb(x,Wa);            
else % Bx <= 150                
if By <= C(3), [Yg rind] = annul(x,d); % C1,C5 < By < C3                
else, [Yg rind] = dispr(x); % C1,C3,C5 < By                
end            
end        
end    
end 
end  

function [Yg rind] = strat(x,Wa) 
rind = 1; Yg = (15400*x/(Wa^0.8))^2;  
end  

function [Yg rind] = wave(x,Wa) 
rind = 2; Yg = 0; 
end  
function [Yg rind] = plug(x,Wa) 
rind = 3; Yg = (27.315*x^0.855 / Wa^0.17)^2; 
end  

function [Yg rind] = slug(x,Wa) 
rind = 4; Yg = (1190*x^0.815/sqrt(Wa))^2; 
end  

function [Yg rind] = bubb(x,Wa) 
rind = 5; Yg = (14.2*x^0.75 / Wa^0.1)^2; 
end  

function [Yg rind] = annul(x,d) 
rind = 6; dx = d; if d > 12, dx = 10; end; 
Yg = ((4.8 - 0.3125*dx)*x^(0.343-0.021*dx))^2; 
end  

function [Yg rind] = dispr(x) 
rind = 7; a0 = 1.4659; a1 = 0.49138; a2 = 0.04887; a3 = -0.000349; 
Yg = (exp(a0 + a1*log(x) + a2*(log(x))^2 + a3*(log(x))^3))^2; 
end 
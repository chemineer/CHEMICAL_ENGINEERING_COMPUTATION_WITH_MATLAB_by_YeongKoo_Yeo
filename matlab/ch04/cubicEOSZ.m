function [Z V] = cubicEOSZ(state,eos,T,P,Tc,Pc,w) 
% Estimation of Z and V using cubic equations of state 
% input 
% state: fluid state (liquid: L, vapor: V) 
% eos: type of equation being used (VDW, RK, SRK, PR) 
% T,P: temperature (K) and pressure (bar) 
% Tc,Pc: critical temperature (K) and pressure (bar) 
% w- acentric factor 
% Tr and Pr (reduced T and P) 
Tr = T/Tc; Pr = P/Pc; nc = length(w); R = 83.14; % cm^3*bar/mol/K 
eos = upper(eos); 
switch eos    
case {'VDW'}         
sm = 0; ep = 0; om = 0.125; ps = 0.42188; mx = [0];    
case{'RK'}         
ep = 0; sm = 1; om = 0.08664; ps = 0.42748;         
al = 1./sqrt(Tr); mx = (Tr.^(-1/4) - 1)./(1 - sqrt(Tr));    
case{'SRK'}         
ep = 0; sm = 1; om = 0.08664; ps = 0.42748;        
al = (1+(0.48+1.574*w-0.176*w.^2).*(1-sqrt(Tr))).^2;         
mc = [0.48 1.574 0.176]; mx = [ones(nc,1) w -w.^2]*mc';    
case{'PR'}         
ep = 1 - sqrt(2); sm = 1 + sqrt(2); om = 0.07780; ps = 0.45724;          
mc = [0.37464 1.54226 0.26992]; mx = [ones(nc,1) w -w.^2]*mc'; 
end 
% calculation of alpha, beta and q 
al = (1 + mx.*(1-sqrt(Tr))).^2; beta = om*Pr./Tr; q = ps*al./(om*Tr); 
% calculation of Z and V 
state = upper(state); 
c(1) = 1; c(2) = (sm +ep)*beta - (1+beta); c(3) = beta*(q + ep*sm*beta -(1+beta)*(sm+ep)); 
c(4) = -beta^2*(q +(1+beta)*ep*sm); Z = roots(c); 
switch state    
case 'V', Z = max(Z);    
case 'L', Z = min(Z); 
end 
V = Z*R*T/P; 
end
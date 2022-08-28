function dPt = dpcatbed(W,Bd,Bl,Pd,Pl,ep,mu,rho) 
% Calculates the pressure drop for flow through a packed bed 
gc = 4.17e8; A = pi*Bd^2/4; G = W/A; % superficial mass flow rate (lb/h/ft^2) 
Ap = (pi*Pd^2/2 + pi*Pd*Pl)/144; Vp = pi*Pl*Pd^2/(4*1728); 
S = Ap*(1-ep)/Vp; ePd = 6*(1-ep)/S; Nre = G*ePd/(2.419*mu*(1 - ep)); 
if Nre <1, fP = 150/Nre; 
elseif Nre < 1e4, fP = 150/Nre + 1.75; 
else, fP = 1.75; end 
dPt = Bl*(1-ep)*G^2*(150*2.419*mu*(1-ep)/ePd/G + 1.75)/(144*ep^3*ePd* gc*rho); 
end 
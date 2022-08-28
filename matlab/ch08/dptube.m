function DPtube = dptube(D,L,rf,v,rho,Visc,Phi,Npass)
% Calculate tube-side pressure drop
% input:
%  D: tube inside diameter (mm)      L: tube length (m)
%  rf: tube roughness (mm)           v- tube-side flow rate (m/s)
%  rho: fluid density (kg/m^3)       Visc: fluid viscosity (Ns/m^2)
%  Phi: viscosity correction factor  Npass: number of tube passes in bundle

% output:
%  DPtube: tube-side pressure drop (N/m^2)
g = 9.81; rfD = rf/D; Dm = 1e-3*D; Nre = Dm*v*rho/Visc;
if Nre <= 2100 % Laminar flow
f = 16/Nre;
elseif Nre <= 4000 % Zigrang & Sylvester friction factor correlation eqn.
t1 = rfD/3.7; t2 = 5.02/Nre; ftm = log10(t1 - t2*log10(t1 + 13/Nre));
f = 1/(4*log10(t1) + t2*ftm)^2;
else % Round friction factor correlation eqn.
f = 1/(3.6*log10(Nre/(0.135*(Nre*rfD + 6.5))))^2;
end
pD = 2*f*rho*v^2*L/Dm; dp1 = pD/Phi; dp2 = 4*rho*L*v^2/(2*g); DPtube = (dp1 + dp2)*Npass;
end

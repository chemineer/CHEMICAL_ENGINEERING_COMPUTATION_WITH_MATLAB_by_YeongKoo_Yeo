function C = rfobj(rf,opdat)
% Retrieve data
we = opdat.we; eta = opdat.eta; S = opdat.S; K = opdat.K;
rhoG = opdat.rhoG; rhoL = opdat.rhoL; rhos = opdat.rhos;
lamb = opdat.lamb; lambs = opdat.lambs; Css = opdat.Css; Cst = opdat.Cst;
F = opdat.F; T = opdat.T; P = opdat.P; alpa = opdat.alpa;
xB = opdat.xB; xD = opdat.xD; xF = opdat.xF;
% Set parameters
rLV = rf/(rf + 1); rDV = 1/(1 + rf);
rpLV = (rf + ((xD-xB)/(xF-xB)))/(1 + rf);
rBV = (((xD-xB)/(xF-xB)) - 1)/(1 + rf);
% Determine the number of trays
N = 1; ye(1) = xD; xe(1) = ye(1)/(alpa*(1-ye(1)) + ye(1));
while xe(N) > xB
if xe(N) > xF
ye(N+1) = rLV*xe(N) + rDV*xD; % rectifying section line
else
ye(N+1) = rpLV*xe(N) - rBV*xB; % stripping section line
end
xe(N+1) = ye(N+1)/(alpa*(1-ye(N+1)) + ye(N+1));
N = N + 1;
end
V = F*(1 + rf)*(xF-xB)/(xD-xB); % vapor rate
d = sqrt((4*V*22.4*760*(T+273.15))/(273*pi*P*K*sqrt((rhoL-rhoG)/rhoL)));
h = 0.6*((N-1)/eta + 1) + 2; % height
w = 14.7*(P/760)*d/2/(we*S - 0.6*14.7*(P/760)) + 0.0032; % steel width
A = 4*pi*(d/2)^2 + pi*d*h; % area
% Cost and objective function
Ceng = lamb*F*(1 + rf)*Css*(xF-xB)/(xD-xB)/lambs; % Energy cost
Ccol = (A + pi*N*(d/2)^2)*w*rhos*Cst; % column cost
C = Ccol/3 + Ceng;
end

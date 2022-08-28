function Htube = htctube(Nre,Pr,D,L,Xk,Phi)
% Calculate convective heat transfer coefficient within tube
% input:
%  Nre: Reynolds number         Pr: Prandtl number
%  D: inside diameter of tube (mm)      L: tube length (m)
%  Xk: thermal conductivity of fluid (W/m/K)    Phi: viscosity correction factor

% output:
%  Htube: convective heat transfer coefficient within tube
Dm = D*1e-3; % mm->m
if Nre <= 2100
Gw = Nre*Pr*Dm/L;
if Gw > 100, Nu = 1.86*Phi*Gw^0.333;
else, Nu = 3.66 + 0.085*Gw*Phi/(1 + 0.047*Gw^0.6667); end
elseif Nre < 1e4, Nu = 0.116*(Nre^0.6667 - 125)*Pr^0.333 * (1+(Dm/L)^0.6667)*Phi;
else, Nu = 0.023*Phi*(Nre^0.8)*(Pr^0.333);
end
Htube = Nu*Xk/Dm;
end



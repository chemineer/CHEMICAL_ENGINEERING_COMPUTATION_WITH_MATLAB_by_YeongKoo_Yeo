function [Z H phi] = phiHeu(x,P,T,state,eos,opdat,mxdat)
% Calculation of mixture enthalpy and fugacity coefficient of component
% Inputs:
%  x: mole fraction vector
%  P, T: pressure (Psia) and temperature (F)
%  state: fluid state ('L': liquid, 'V': vapor)
%  eos: equation of state ('RK', 'SRK', or 'PR')
%  mxp: physical property structure
% mxp.Pc, mxp.Tc: critical pressure (Psia) and temperature (R) vector
% mxp.w: acentric factor for each component (vector)
% mxp.k: nxn diagonal matrix of binary interaction parameters
% mxp.Afi: heat capacity integration coefficients for ideal gases
% Outputs:
%  Z: compressibility factor
%  H: mixture enthalpy (Btu/lbmol)
%  phi: fugacity coefficient
w = mxdat.w; k = mxdat.k; Tc = mxdat.Tc/1.8; Pc = 6894.8*mxdat.Pc; % T and P in K and Pa
T = (T-32)/1.8 + 273.15; P = 6894.8*P; Afi = mxdat.Afi; % Btu/lbmole(1Btu/lbmole =2.326 J/mol)
nc = opdat.nc; R = 8.314; % gas constant: m^3 Pa/(mol K) = J/mol-K
Tr = T./Tc; Pr = P./Pc; nc = opdat.nc; eos = upper(eos); state = upper(state);
switch eos
case{'RK'}
ai = sqrt(0.4278./(Pc.*Tr.^2.5)); bi = 0.0867./(Pc.*Tr); A = sum(x.*ai); B = sum(x.*bi);
Z = roots([1 -1 B*P*(A^2/B-B*P-1) -A^2*(B*P)^2/B]);
case{'SRK'}
mx = 0.48+1.574*w-0.176*w.^2; al = (1+mx.*(1-sqrt(Tr))).^2;
ai = 0.42747*al.*Pr./(Tr.^2); bi = 0.08664*Pr./Tr;
am = sqrt(ai'*ai).*(1-k); A = x*am*x'; B = sum(x.*bi); Z = roots([1 -1 A-B-B^2 -A*B]);
case{'PR'}
mx = 0.37464+1.54226*w-0.26992*w.^2; al = (1+mx.*(1-sqrt(Tr))).^2;
ai = 0.45723553*al.*Pr./(Tr.^2); bi = 0.0777961*Pr./Tr;
am = sqrt(ai'*ai).*(1-k); A = x*am*x'; B = sum(x.*bi);
Z = roots([1 B-1 A-3*B^2-2*B B^3+B^2-A*B]);
end
iz = abs(imag(Z)); Z(and(iz>0,iz<=1e-6)) = real(Z(and(iz>0,iz<=1e-6)));
for i = 1:length(Z), zind(i) = isreal(Z(i)); end
Z = Z(zind);
if state == 'L', Z = min(Z);
else, Z = max(Z); end
V = R*T*Z/P; % m^3/mol
% Enthalpy and fugacity coefficients
% Hv0 = int(Cp) Btu/lbmol, Afi uses T in F
Te = (T-273.15)*1.8 + 32; % T: K->F
Hv0 = x*(Afi(:,1)*Te + Afi(:,2)*Te^2/2 + Afi(:,3)*Te^3/3 +Afi(:,4)*Te^4/4 + Afi(:,5)*Te^5/5);
Hv0 = Hv0*2.326; % 1Btu/lbmole=2.326J/mol
switch eos
case{'RK'}
H = Hv0 + R*T*(Z - 1 - 3*(A^2)*log(1 + B*P/Z)/(2*B));
phi = exp((Z-1)*bi/B - log(Z-B*P) - (A^2/B)*(2*ai/A - bi/B)*log(1 + B*P/Z));
case{'SRK'}
hsum = 0;
for i = 1:nc
for j = 1:nc
hsum = hsum + x(i)*x(j)*am(i,j)*(1 - mx(i)*sqrt(Tr(i))/(2*sqrt(al(i))) -...
mx(j)*sqrt(Tr(j))/(2*sqrt(al(j))));
end
end
H = Hv0 + R*T*(Z - 1 - log((Z + B)/Z)*hsum/B);
phi = exp((Z-1)*bi/B - log(Z-B) - (A/B)*(2*sqrt(ai)/sqrt(A) - bi/B)*log((Z+B)/Z));
case{'PR'}
hsum = 0;
for i = 1:nc
for j = 1:nc
hsum = hsum + x(i)*x(j)*am(i,j)*(1 - mx(i)*sqrt(Tr(i))/ (2*sqrt(al(i))) -...
mx(j)*sqrt(Tr(j))/(2*sqrt(al(j))));
end
end
H = Hv0 + R*T*(Z - 1 - log((Z + B)/Z)*hsum/B);
phi = exp((Z-1)*bi/B - log(Z-B) - (A/B/sqrt(8))*(2*x(1:end)*am(1:end,:)/A -...
bi/B)*log((Z+(1+sqrt(2))*B)/(Z+(1-sqrt(2))*B)))
end
H = H/2.326; % J/mol -> Btu/lbmol
end

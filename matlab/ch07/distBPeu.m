function [x,y,T,L,V,iter] = distBPeu(opdat,mxdat)
% Calculation of multicomponent distillation using BP method (English unit)
% Set data
eos = opdat.eos; N = opdat.N; nc = opdat.nc; F = opdat.F; Tf = opdat.Tf; Pf = opdat.Pf; hF = opdat.hF;
P = opdat.P; T = opdat.T; V = opdat.V; L = opdat.L; criv = opdat.criv;
U = opdat.U; W = opdat.W; Q = opdat.Q; T0 = opdat.T0; nf = opdat.nf;
hV = opdat.hV; hL = opdat.hL; fstate = upper(opdat.fstate);
x = opdat.x; y = opdat.y; z = opdat.z; Pc = mxdat.Pc; Ant = mxdat.Ant;
% Initialization
a = zeros(1,N); b = zeros(1,N); d = zeros(1,N); V(2) = L(1) + V(1) + U(1) + W(1) - F(1);
% Feed enthalpy
for j = 1:length(nf)
[Z H phi] = phiHeu(z(:,nf(j))',Pf(nf(j)),Tf(nf(j)),fstate,eos,opdat, mxdat); hF(nf(j)) = H; % Btu/mol
end
T = T0; Told = T; criT = 10; iter = 1;
% Initialization of K(i,j): K(i,j) at T0
for j = 1:N
Pv0(:,j) = Pc(:).*exp(Ant(:,1) - Ant(:,2)./(T(j) + Ant(:,3))); K(:,j) = Pv0(:,j)./P(j); % Raoult's law
end
while criT > criv
% Calculate x(i,j)
for i = 1:nc
for j = 2:N-1
Pj(j) = V(j) - V(1) + sum(F(1:j-1) - U(1:j-1) - W(1:j-1));
Qj(j) = V(1) - (V(j) + W(j))*K(i,j) - U(j) - V(j+1) - sum(F(1:j) - U(1:j) - W(1:j));
Rj(j) = V(j+1)*K(i,j+1); Sj(j) = -F(j)*z(i,j);
end
Pj(N) = V(N) - V(1) + sum(F(1:N-1) - U(1:N-1) - W(1:N-1));
Qj(1) = V(1) - (V(1) + W(1))*K(i,1) - U(1) - V(2) - (F(1) - U(1) - W(1));
Qj(N) = V(1) - (V(N) + W(N))*K(i,N) - U(N) - sum(F(1:N) - U(1:N) - W(1:N));
Rj(1) = V(2)*K(i,2); Sj(1) = -F(1)*z(i,1); Sj(N) = -F(N)*z(i,N);
r(1) = Rj(1)/Qj(1); s(1) = Sj(1)/Qj(1);
for j = 2:N-1
r(j) = Rj(j)/(Qj(j) - r(j-1)*Pj(j)); s(j) = (Sj(j) - s(j-1)*Pj(j))/(Qj(j) - r(j-1)*Pj(j));
end
s(N) = (Sj(N) - s(N-1)*Pj(N))/(Qj(N) - r(N-1)*Pj(N)); x(i,N) = s(N);
for j = N-1:-1:1, x(i,j) = s(j) - r(j)*x(i,j+1); end
end
for j = 1:N, x(:,j) = x(:,j)/sum(x(:,j)); end
% Calculate new T(j) and y(i,j) using BP method
for j = 1:N % calculate y(i,j)
y(:,j) = K(:,j).*x(:,j); y(:,j) = y(:,j)/sum(y(:,j));
end
% Antoine eqn. (T: deg.F)
for j = 1:N % calculate T(j) (T in F)
f = @(Tv) sum(Pc'.*exp(Ant(:,1) - Ant(:,2)./(Tv +Ant(:,3))).*x(:,j)/P(j)) - 1; T(j) = fzero(f,T(j));
end
for j = 1:N
lstate = 'L';
[Z H phi] = phiHeu(x(:,j)',P(j),T(j),lstate,eos,opdat,mxdat);
phiL = phi; hL(j) = H; % Btu/lbmol
vstate = 'V'; [Z H phi] = phiHeu(y(:,j)',P(j),T(j),vstate,eos,opdat, mxdat);
phiV = phi; hV(j) = H; % Btu/lbmol
K(:,j) = phiL./phiV;
end
% Calculate loads for the condenser and reboiler
Q(1) = (L(1) + V(1) + U(1) + W(1) - F(1))*hV(2) - V(1)*hV(1) - (L(1) + U(1))*hL(1) + F(1)*hF(1);
Q(N) = sum(F.*hF - U.*hL - W.*hV) - sum(Q(1:N-1)) - V(1)*hV(1) - L(N)*hL(N);
% Calculate V(j) and L(j) sequentially
for j = 2:N
a(j) = hL(j-1) - hV(j); b(j-1) = hV(j) - hL(j-1);
d(j) = (hL(j)-hL(j-1))*sum(F(1:j-1)-U(1:j-1)-W(1:j-1)) + ...
F(j)*(hL(j)-hF(j)) + W(j)* (hV(j)-hL(j)) + Q(j);
end
a(1) = -hV(1); b(N) = -hL(N); d(1) = F(1)*(hL(1)-hF(1)) + W(1)*(hV(1)-hL(1)) + Q(1);
for j = 2:N-1
V(j+1) = (d(j) - a(j)*V(j))/b(j); L(j) = V(j+1) - V(1) + sum(F(1:j) - U(1:j) - W(1:j));
end
% Check convergence
criT = sum(abs(Told - T)); Told = T; iter = iter + 1;
end
end

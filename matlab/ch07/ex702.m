% mdif.m
Ca1=2.229e-4; Cb1=0; Cc1=7.208e-3; Ca2=0; Cb2=2.701e-3; Cc2=4.73e-3;
D1=1.075e-4; D2=1.245e-4; D3=1.47e-4; P=0.2;T=328; R=82.057e-3; Ct = P/(R*T);
L=0.001; Nb=-4.143e-4; Nc=0; Na=-D3*(Ca2-Ca1)/L;
zspan = [0 L]; c0 = [Ca1 Cb1 Cc1]; critN = 1e-10; errA = 1; Na1 = Na/2; Na2 = 2*Na;
iter = 1;
while errA > critN
Nam = (Na1+Na2)/2;
[z c1] = ode45(@mf,zspan,[Ca1,Cb1,Cc1],[],D1,D2,D3,Na1,Nb,Nc,Ct);
[z c2] = ode45(@mf,zspan,[Ca1,Cb1,Cc1],[],D1,D2,D3,Na2,Nb,Nc,Ct);
[z cm] = ode45(@mf,zspan,[Ca1,Cb1,Cc1],[],D1,D2,D3,Nam,Nb,Nc,Ct);
if c1(end,1)*cm(end,1) < 0, Na2 = Nam; % check whether Ca(z2)=0 is satisfied
else, Na1 = Nam; end
errA = abs(Na1 - Na2); iter = iter+1;
end
c = cm; xa = c(:,1)/Ct; xb = c(:,2)/Ct; xc = c(:,3)/Ct;
plot(z,xa,z,xb,':',z,xc,'.-'), legend('x_A','x_B','x_C'), xlabel('Distance z(m)'), ylabel('Mole fraction')
iter, Na = Nam

function dc = mf(z,c,D1,D2,D3,Na,Nb,Nc,Ct)
% c1=Ca, c2=Cb, c3=Cc, D1=Dab, D2=Dbc, D3=Dac
xa = c(1)/Ct; xb = c(2)/Ct; xc = c(3)/Ct;
dc = [(xa*Nb-xb*Na)/D1 + (xa*Nc-xc*Na)/D3;
(xb*Na-xa*Nb)/D1 + (xb*Nc-xc*Nb)/D2;
(xc*Na-xa*Nc)/D3 + (xc*Nb-xb*Nc)/D2];
end

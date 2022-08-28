% rdist.m
% s1: T0, s2: T1, s3: Tf, s4: T2, s5: T3, s6: x11, s7: x21
% s8: x12, s9: x22, s10: x13, s11: x23, s12: V1, s13: V2, s14: V3
F = 1; D = 0.25; z = [0.23 0.77]; P = 760*120/14.7; Q = 1e4;
T0 = 200; T1 = 145; Tf = 200; T2 = 190; T3 = 210; x11 = 0.65; x21 = 0.35;
x12 = 0.43; x22 = 0.57; x13 = 0.33; x23 = 0.76; V1 = 1.1; V2 = 1; V3 = 1.1;
s0 = [T0 T1 Tf T2 T3 x11 x21 x12 x22 x13 x23 V1 V2 V3];
s = fsolve(@fd,s0,[],F,D,z,P,Q);
T0 = s(1); T1 = s(2); Tf = s(3); T2 = s(4); T3 = s(5);
x11 = s(6); x21 = s(7); x12 = s(8); x22 = s(9); x13 = s(10); x23 = s(11); V1 = s(12);
V2 = s(13); V3 = s(14);
fprintf('T0=%8.4f, T1=%8.4f, Tf=%8.4f, T2=%8.4f, T3=%8.4f \n',T0,T1,Tf,T2,T3);
fprintf('x11=%6.4f, x12=%6.4f, x13=%6.4f\n', x11, x12, x13);
fprintf('x21=%6.4f, x22=%6.4f, x23=%6.4f\n', x21, x22, x23);
fprintf('V1=%6.4f, V2=%6.4f, V3=%6.4f\n', V1, V2, V3);
function f = fd(s,F,D,z,P,Q)
% i=1: n-butane, i=2:n-Pentane)
% s1: T0, s2: T1, s3: Tf, s4: T2, s5: T3, s6: x11, s7: x21
% s8: x12, s9: x22, s10: x13, s11: x23, s12: V1, s13: V2, s14: V3
% Antoine coefficients and parameters
A = [6.80776 6.85296]; B = [935.77 1064.84]; C = [238.789 232.012];
h1c = [0.04 29.6]; h2c = [0.025 38.5]; H1c = [-0.04 43.8 8003]; H2c = [0.007 31.7 12004];
% Equilibrium constants
k0 = 10.^(A - B./((s(1)-32)*5/9 + C)) / P; k1 = 10.^(A - B./((s(2)-32)*5/9 + C)) / P;
kf = 10.^(A - B./((s(3)-32)*5/9 + C)) / P; k2 = 10.^(A - B./((s(4)-32)*5/9 + C)) / P;
k3 = 10.^(A - B./((s(5)-32)*5/9 + C)) / P;
% Composition vector
x1 = [s(6) s(7)]; x2 = [s(8) s(9)]; x3 = [s(10) s(11)]; x0 = k1.*x1;
% Enthalpy
hL0 = sum(([s(1)^2 s(1)]*[h1c' h2c']).*x0); hL1 = sum(([s(2)^2 s(2)]*[h1c' h2c']).*x1);
hLf = sum(([s(3)^2 s(3)]*[h1c' h2c']).*z); hL2 = sum(([s(4)^2 s(4)]*[h1c' h2c']).*x2);
hL3 = sum(([s(5)^2 s(5)]*[h1c' h2c']).*x3); hV1 = sum(([s(2)^2 s(2) 1]*[H1c' H2c']).*k1.*x1);
hV2 = sum(([s(4)^2 s(4) 1]*[H1c' H2c']).*k2.*x2); hV3 = sum(([s(5)^2 s(5) 1]*[H1c' H2c']).*k3.*x3);
% Mass balance
B = F - D; L0 = s(12) - D; L1 = s(13) - D; L2 = s(14) + B; L3 = B;
% Definition of equations
f(1,1) = sum(k0.*x0) - 1; f(2,1) = sum(k1.*x1) - 1; f(3,1) = sum(kf.*z) - 1; f(4,1) = sum(k2.*x2) - 1;
f(5,1) = sum(k3.*x3) - 1; f(6,1) = -((s(12) - L0)*k1(1) + L1)*s(6) + s(13)*k2(1)*s(8);
f(7,1) = -((s(12) - L0)*k1(2) + L1)*s(7) + s(13)*k2(2)*s(9);
f(8,1) = -s(12)*hV1 + s(13)*hV2 - L1*hL1 + L0*hL0;
f(9,1) = L1*s(6) - (s(13)*k2(1)+L2)*s(8) + s(14)*k3(1)*s(10) + F*z(1);
f(10,1) = L1*s(7) - (s(13)*k2(2)+L2)*s(9) + s(14)*k3(2)*s(11) + F*z(2);
f(11,1) = -s(13)*hV2 + s(14)*hV3 + hLf + L1*hL1 - L2*hL2;
f(12,1) = L2*s(8) - (s(14)*k3(1) + B)*s(10);
f(13,1) = L2*s(9) - (s(14)*k3(2) + B)*s(11);
f(14,1) = -s(14)*hV3 + Q + L2*hL2 - L3*hL3;
end

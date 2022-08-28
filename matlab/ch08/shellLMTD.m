% shellLMTD.m : number of required shells and LMTD
clear all;
% Data
T1 = input('Hot fluid inlet temperature (deg.F): ');
T2 = input('Hot fluid outlet temperature (deg.F): ');
t1 = input('Cold fluid inlet temperature (deg.F): ');
t2 = input('Cold fluid outlet temperature (deg.F): ');
N = 1; Dt1 = T2 - t1; Dt2 = T1 - t2;
P = (t2 - t1)/(T1 - t1); R = (T1 - T2)/(t2 - t1); % Compute P and R
A = -1; B = -1; F = 0.1; % Initial F
while (F <= 0.75) % Test F and R
if ( R > 1 || R < 1)
while (A < 0)
Pp = (1 - ((P*R-1)/(P-1))^(1/N)) / (R-((P*R-1)/(P-1))^(1/N));
A = (2/Pp - 1 - R + sqrt(R^2 +1)) / (2/Pp - 1 - R - sqrt(R^2 +1));
if (A < 0), N = N + 1; end
end
F = ( sqrt(R^2 +1) * log10((1-Pp)/(1-Pp*R)) / (R-1) ) / (log10(A));
else  % R = 1
while (B < 0)
Ppp = P / (N-P*(N-1)); B = (2/Ppp - 2 + sqrt(2)) / (2/Ppp - 2 - sqrt(2));
if (B < 0), N = N + 1; end
end
F = ( sqrt(R^2 +1) * Ppp / (log(10*(1-Ppp))) ) / (log10(B));
end
if (F <= 0.75), N = N + 1; end
end
if (Dt1 == Dt2), LMTD = Dt1; else, LMTD = (Dt1-Dt2)/(log(Dt1/Dt2)); end %Compute LMTD
cLMTD = F*LMTD; % Corrected LMTD
% Output
fprintf('\nNumber of shells = %3d\n', N); fprintf('F factor = %9.4f\n', F);
fprintf('Corrected LMTD = %9.4f\n', cLMTD);

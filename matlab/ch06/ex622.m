% pbrpar.m 
r = 1e-10*[71.0 71.3 41.6 19.7 42.0 17.1 71.8 142.0 284.0 47.0 71.3 117.0 127.0 131.0 133.0 41.8]; 
Pt = [1 1 1 1 1 1 1 1 1 0.5 1 5 10 15 20 1]; Ph = [1 1 1 1 1 1 1 2 4 1 1 1 1 1 1 1]; 
Pm = [1 4 0 0 1 0 0 0 0 0 0 0 0 0 0 1]; Pb = [0 0 1 4 1 5 0 0 0 0 0 0 0 0 0 1]; 
n = length(r); A = [ones(n,1) Pb' Pt']; b = Ph.*Pt./r; 
x = inv(A'*A)*A'*b'; k = 1/x(1), Kb = x(2)>k, Kt = x(3)*k 
rc = k*Ph.*Pt./(1 + Kb*Pb + Kt*Pt); nc = 1:n; 
plot(nc,r*1e10,'o',nc,rc*1e10,'*'), legend('Data','Estimated'), xlabel('Run #'), ylabel('-10^{10}r') 
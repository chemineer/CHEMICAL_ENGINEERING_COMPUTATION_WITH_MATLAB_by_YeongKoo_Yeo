% flashdrum.m: flash calculation using modified Raoult's law 
zf = [0.1 0.25 0.5 0.15];  % feed composition 
A = [6.64380 6.82915 6.80338 6.80776]; % Antoine constants 
B = [395.74 663.72 804.00 935.77]; 
C = [266.681 256.681 247.040 238.789]; 
P = [16 18 20 24]*760; % pressure (mmHg) 
% compute alpha 
T = 50; Pv = 10.^(A - B./(C + T)); % vapor pressure (mmHg) (T: deg.C) 
x = []; y = []; alpha = []; n = length(P); 
for i = 1:n    
k = Pv/P(i); ax = fzero(@sumf,0.5,[],zf,k); xv = zf./(1 + ax*(k - 1));    
yv = k.*xv; x = [x xv']; y = [y yv']; alpha = [alpha ax]; 
end 
Patm = P/760; 
% dew-point: alpha = 1 
Tdew = []; 
for i = 1:n, T = fzero(@zdk,50,[],zf,P(i),A,B,C); Tdew = [Tdew T]; end 
% bubble-point: alpha = 0 
Tbub = []; 
for i = 1:n, T = fzero(@zmk,50,[],zf,P(i),A,B,C); Tbub = [Tbub T]; end 
Patm, x, y, alpha, Tdew, Tbub  

function f = sumf(alp,zf,k) 
f = sum(zf.*(1-k)./(1 + alp*(k-1))); 
end 
function f = zdk(T,zf,P,A,B,C) % dew-point calculation: alpha=1 
k = 10.^(A - B./(C + T)) / P; f = sum(zf./k) - 1; 
end 
function f = zmk(T,zf,P,A,B,C) % bubble-point calculation: alpha=0 
k = 10.^(A - B./(C + T)) / P; f = sum(zf.*k) - 1; 
end 
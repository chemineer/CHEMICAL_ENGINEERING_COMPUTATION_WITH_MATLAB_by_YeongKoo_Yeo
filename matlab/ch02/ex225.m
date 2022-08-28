% antvp.m: vapor pressure by Antoine eqn. 
T =  [25 27 30 31 35 36 37]; Pv =  [15.8 26.43 39.56 94.91 428.35 861.71 1851.24]; 
X  =  [ones(1,length(T)); 1./T; log10(Pv)./T]'; 
b  =  inv(X'*X)*X'*log10(Pv)'; A  =  b(1), C  =  -b(3), B  =  b(2)-A*C 
Ti  =  T(1):0.1:T(end); Pvi  =  10.^(A + B./(Ti+C)); 
plot(Ti,Pvi,T,Pv,'o'), xlabel('T(C)'), ylabel('Pv(mmHg)'), legend('Fitting','Data') 
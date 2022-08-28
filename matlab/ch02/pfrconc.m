function dC  =  pfrconc(t,C,pf) 
% Retrieve data 
k  =  pf.k; v  =  pf.v; C0  =  pf.C0; L  =  pf.L; n  =  pf.n; 
% Initialization 
h  =  L/n; dC  =  zeros(n,1); 
% Difference equations 
for m  =  1:n 
if m  ==  1, s  =  (v/2/h)*(C(m+1) - C0); 
elseif m == n, s  =  (v/h)*(C(m) - C(m-1)); 
else, s  =  (v/2/h)*(C(m+1) - C(m-1)); end 
dC(m)  =  -s - k*C(m); 
end 
end 
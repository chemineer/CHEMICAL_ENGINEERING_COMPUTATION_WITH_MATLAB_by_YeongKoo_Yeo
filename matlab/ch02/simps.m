function z  =  simps(f,a,b,n) 
% Implements Simpson 1/3 rule to integrate f(x) from a to b 
% input 
% f: function handle to be integrated 
% a,b: limits of integration 
% n: number of subinterval points (x1,??,xn) 
% Output 
% z: integral of f(x) 
if ceil(n/2) - floor(n/2) < 1, n  =  n+1; end % n must be odd 
h  =  (b-a)/(n-1); s  =  f(a); 
if n <=  2, display('Too few subintervals.'); return; end 
if n  ==  3, z  =  h*(s + 4*f((a+b)/2) + f(b))/3; return; end 
for k  =  2:2:n-1, x  =  a + h*(k-1); s  =  s + 4*f(x); end 
for k  =  3:2:n-2, x  =  a + h*(k-1); s  =  s + 2*f(x); end 
s  =  s + f(b); z  =  h*s/3; 
end 
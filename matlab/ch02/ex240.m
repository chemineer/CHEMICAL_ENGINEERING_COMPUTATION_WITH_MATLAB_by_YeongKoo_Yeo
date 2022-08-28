% filt.m 
t  =  [4.4 9.5 16.3 24.6 34.7 46.1 59.0 73.6 89.4 107.3]; % t (sec) 
V  =  1e-3*[0.498 1.0 1.501 2.0 2.498 3.002 3.506 4.004 4.502 5.009]; % V (m^3) 
A  =  0.04; cs  =  20; vis  =  8.937e-4; dp  =  3e5; % data 
dtV  =  diff(t)./diff(V); % numerical differentiation using diff 
n  =  length(V); Vm  =  (V(1:n-1) + V(2:n))/2; % midpoint of each interval 
k1  =  vis*cs/(A^2*dp); k2  =  vis/(A*dp); % define k1 and k2 
Y  =  dtV'; X  =  [k1*Vm; k2*ones(1,n-1)]'; C  =  inv(X'*X)*X'*Y; % regression 
fprintf('alpha  =  %g, Rm  =  %g\n', C(1), C(2)) 
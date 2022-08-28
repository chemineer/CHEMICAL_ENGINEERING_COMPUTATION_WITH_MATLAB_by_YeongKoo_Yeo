% numint.m 
t  =  [20.1 49.8 81.5 109.8 140.2 171.1 201.4 229.5]; % temp.(deg.C) 
Cp  =  [28.98 29.11 29.96 29.47 29.67 29.91 30.02 30.14]; % Cp (J/mol/deg.C) 
n  =  6.5; t1  =  55; t2  =  185; % data 
% step 1: find Cp(t1) and Cp(t2) by interpolation using piecewise cubic spline method. 
Cp1 = interp1(t, Cp, t1,'spline'); Cp2 = interp1(t, Cp, t2, 'spline'); % Cp at t =  t1 and t  =  t2 
% step 2: numerical interpolation by trapz on the subintervals [t1 t2] and [Cp1 Cp2] 
ts  =  [t1 81.5 109.8 140.2 171.1 t2]; % subinterval of t: t1 < =  t < =  t2 
Cps  =  [Cp1 29.96 29.47 29.67 29.91 Cp2]; % subinterval of Cp: Cp1 < =  Cp < =  Cp2 
delH  =  n*trapz(ts,Cps) % perform numerical integration 
function dfv = dfun(x)
df = [2 2*x(2)]; % gradient of objective function
dh = [2*x(1) 2*x(2)]; % gradient of equality constraint (h(x) = 0)
dg = [1 0;-1 0;0 1;0 -1]; % gradient of inequality constraint (gi(x) >= 0)
dfv = [df' dh' dg'];
end

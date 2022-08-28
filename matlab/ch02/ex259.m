% Motion of a Vibrating String
f1 = @(x) x.*(1-x); f2 = @(x) 0; g0 = @(t) 0; g1 = @(t) 0; xspan = [0 1]; tspan =  [0 1]; 
nx  =  20; nt  =  40; alpa  =  1; [u r]  =  hypbpde(f1,f2,g0,g1,xspan,tspan,nx, nt,alpa);

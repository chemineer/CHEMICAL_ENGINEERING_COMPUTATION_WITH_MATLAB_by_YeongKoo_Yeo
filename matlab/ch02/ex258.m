% Temperature Distribution in a Rod
f = @(x) x.^3 - 2*x.^2 + 1.5*x; g0 = @(t) 0; g1 = @(t) 2; nx = 10; nt = 50; alpa =  0.8; tf  =  0.1; 
u  =  parapde(f,g0,g1,tf,nx,nt,alpa);
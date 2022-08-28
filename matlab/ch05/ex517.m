% Flow in a Pipeline Network
L = 100*[1 3 12 3 12 12 3]; dP0 = -15e5; rho = 997.08; mu = 8.931e-4; D = 0.154; x0 = 0.1*ones(1,7);  
x = fsolve(@pnetq,x0,[],D,rho,mu,dP0,L)
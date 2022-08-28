% Internal Diameter of a Pipe
W=0.36; rho=85.9; mu=4.4e-4; rf=5e-6; dP=1e4; L=4; K=0.3; D0=0.1;
D = fzero(@dwfun,D0,[],W,rho,mu,rf,dP,L,K)
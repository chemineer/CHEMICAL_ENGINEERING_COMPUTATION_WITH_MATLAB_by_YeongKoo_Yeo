% Pipe Diameter for Non-Newtonian Flow
W=6.67; rho=961; dP=1.5e4; rf=5e-6; L=10; k=1.48; n=0.64; K=1.8; 
D0 = 0.1; D = fzero(@nnDfun,D0,[],W,rho,dP,rf,L,k,n,K) 
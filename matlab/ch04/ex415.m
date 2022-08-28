% actunifgam.m 
Nc = 2; k = 3; % numbers of components(Nc) and functional groups(k) 
R = [0.9011 0.6744 1.2070]; % vector of volumes for each functional group 
Q = [0.8480 0.5400 0.9360]; % vector of surface areas for each functional group 
nu = [2 2; 1 5; 1 0]; % number of functional groups (row: functional groups k, column: components i) 
amn = [0 0 255.7; 0 0 255.7; 65.33 65.33 0]; % matrix of group interaction parameters 
T = 308.15; x = [0.4 0.6]; % mole fraction of each component (T: K) 
gam = unifgam(k,R,Q,nu,amn,Nc,x,T)
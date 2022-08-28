% actmixunifgam.m 
Nc = 4; k = 6; % numbers of components(Nc) and functional groups(k) 
R = [0.9011 0.6744 0.4469 0.2195 0.5313 1.0000]; % vector of volumes for each functional group 
Q = [0.848 0.540 0.228 0.000 0.400 1.200]; % vector of surface areas for each functional group 
% number of functional groups contained in each component 
% (row: functional groups k, column: components i) 
nu = [2 1 3 0; 4 1 1 0; 0 0 1 0; 0 0 1 0; 0 0 0 6; 0 1 0 0];  
amn = [0 0 0 0 61.13 986.5; 0 0 0 0 61.13 986.5; 0 0 0 0 61.13 986.5;...           
0 0 0 0 61.13 986.5; -11.12 -11.12 -11.12 -11.12 0 636.1;...           
156.4 156.4 156.4 156.4 89.60 0]; % matrix of group interaction parameters 
T = 334.82; x = [0.162 0.068 0.656 0.114]; % % mole fraction of each component (T: K) 
gam = unifgam(k,R,Q,nu,amn,Nc,x,T)       
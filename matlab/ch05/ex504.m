% Terminal Velocity of a Falling Particle
rp = 1780; ro = 994.6; dp = 2e-4; mu = 8.931e-4; vt0 = 1e-3; 
vt = fzero(@vtfun,vt0,[], rp, ro, mu, dp) 
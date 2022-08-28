% mflux.m 
D  =  3.26e-8; x  =  -1e-3:1e-4:1e-3; 
C  =  @(x) -1.3e6*x.^4 + 7.1e3*x.^3 - 14*x.^2 - 0.364*x + 0.001; 
f  =  C(x); numf  =  diff(f)./diff(x); % numerical differentiation 
NCO2  =  -D*numf; % flux of CO2 
nfl  =  NCO2(1); nfr  =  NCO2(end); netf  =  nfl - nfr; % net flux 
fprintf('Flux at x = -0.001: %g, Flux at x = 0.001: %g, Net flux: %g\n', nfl, nfr, netf) 
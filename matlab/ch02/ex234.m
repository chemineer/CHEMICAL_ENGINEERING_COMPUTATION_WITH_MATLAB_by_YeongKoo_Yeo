% cp1dintp.m: 1-D interpolation for Cp data 
T  =  373:100:873; Cp  =  [29.189 29.291 29.462 29.678 29.971 30.269]; % data 
lCp  =  interp1(T, Cp, 580, 'linear'); pCp  =  interp1(T, Cp, 580, 'pchip'); 
sCp  =  interp1(T, Cp, 580, 'spline'); nCp  =  interp1(T, Cp, 580, 'nearest'); 
fprintf('linear: %6.4f\npchip: %6.4f\nspline: %6.4f\nnearest: %6.4f\n', lCp,pCp,sCp,nCp);  
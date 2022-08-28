% pbfun.m: normal probability density function 
% 'pkg load statistics' required
clear all;
mu  =  0; sigma  =  1; lx  =  -0.8; ux  =  0.8; zl  =  2; 
fx  =  @(x) normpdf(x,mu,sigma); % normal probability density function f(x) 
P  =  quad (fx,lx,ux); % probability assuming any value between lx and ux 
Pr  =  normcdf(ux,mu,sigma) - normcdf(lx,mu,sigma); % probability of z lying between lx and ux 
Pz  =  1 - normcdf(zl, mu,sigma); % probability of observing z >= zl 
fprintf('Probability of z assuming any value between %g and %g = %g \n',lx,ux,P); 
fprintf('Probability of z lying between %g and %g  =  %g\n',lx,ux,Pr); 
fprintf('Probability of observing z > =  %g  =  %g\n',zl,Pz);  


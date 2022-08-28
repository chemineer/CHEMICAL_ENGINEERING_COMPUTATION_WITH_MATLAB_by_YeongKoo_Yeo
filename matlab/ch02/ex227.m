% diffwood.m: estimation of water diffusivity in the wood 
% 'pkg install -forge optim' required
% 'pkg load optim' required
clear all; t  =  [5:5:30]; % time(hr) 
x  =  [0.245, 0.230, 0.211, 0.197, 0.187, 0.176]; % free moisture in the wood 
z  =  0.03; % wood thickness 
f  =  @(C,t) 8*C(1)*(exp(-C(2)*t*(pi/2/z)^2) + exp(-9*C(2)*t*(pi/2/z)^2)/ 9)/pi^2; 
C0  =  [0.3; 0]; % initial guess (C(1) = x0, C(2) = D) 
C  =  nlinfit(t,x,f,C0); % nonlinear regression by built-in function nlinfit 
x0  =  C(1); D  =  C(2); tv  =  0:0.1:30; 
fprintf('Initial free moisture content of the wood  =  %g\n', x0) 
fprintf('Diffusivity of water in the wood  =  %g\n', D) 
plot(tv,f(C,tv),t,x,'o'), grid, xlabel('t(hr)'), ylabel('x(kg H_2O/ kg wood)') 
legend('Fitting by built-in fun nlinfit','Data','location','best')
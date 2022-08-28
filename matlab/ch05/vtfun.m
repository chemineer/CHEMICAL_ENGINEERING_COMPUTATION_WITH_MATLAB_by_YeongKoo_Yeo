function fvt = vtfun(x, rp, ro, mu, dp) 
% x: unknown terminal velocity 
g = 9.8; Nre = dp*ro*x/mu; % Reynolds number 
if Nre < 0.1, Cd = 24./Nre; 
elseif Nre <1000, Cd = (1 + 0.14*Nre.^0.7)*24./Nre; 
elseif Nre < 3.5e5, Cd = 0.44; 
else, Cd = 0.19 - 8e4./Nre; end 
fvt = 3*Cd*ro*x.^2 - 4*g*(rp - ro)*dp; % nonlinear equation 
end 
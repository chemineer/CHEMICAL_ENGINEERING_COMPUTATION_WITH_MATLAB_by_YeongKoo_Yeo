function v = hzpipe(L,R,mu,delP) 
% Average velocity and velocity distribution of the laminar flow inside a horizontal circular pipe 
% inputs 
% L: length of pipe (m), R: radius of pipe (m) 
% mu: viscosity (kg/m/sec), delP: pressure drop (Pa) 
% output 
% v: average velocity (m/s) 
r = linspace(-R,R,100);  
vr = delP*R^2.*(1 - (r/R).^2)/(4*mu*L); v = delP*R.^2/(8*mu*L);  
plot(vr,r), grid, axis([0 1.5 -R,R]), xlabel('v_x(r)'), ylabel('r'), axis([0 1 -R R]) 
end 
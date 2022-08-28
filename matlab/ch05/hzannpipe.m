function avgv = hzannpipe(L,R1,R2,mu,delP) 
% Calculates the average velocity and plots velocity and shear stress  
% distribution for Newtonian laminar flow in a horizontal annulus. 
% input 
% L: length of pipe (m) 
% R1,R2: inside radius of inner and outer pipe (m) 
% mu: viscosity (kg/m/s),  delP: pressure drop (Pa) 
% output 
% avgv: average velocity (m/s) 
n = 100; r = linspace(R1,R2,n); h = (R2-R1)/n; 
v = delP*(R2^2 -r.^2 + (R2^2 -R1^2)/(log(R2/R1)).*log(r/R2))./(4*mu*L); 
avgv = delP*(R1.^2 + R2.^2 - (R2.^2 -R1.^2)./(log(R2./R1)))./(8*mu*L); 
dv = diff(v); dv = [dv dv(end)]; taurx = -mu.*dv./h; 
subplot(1,2,1), plot(v,r), grid, axis([0 0.2 R1 R2]), xlabel('v_x(r)'), 
ylabel('r') 
subplot(1,2,2), plot(taurx,r), grid, axis tight, xlabel('\tau_{rx}(r)'), 
ylabel('r') 
end 
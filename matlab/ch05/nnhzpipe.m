function [avgv taurx] = nnhzpipe(L,R,K,n,delP) 
% average velocity and velocity profile for a non-Newtonian fluid flowing in a horizontal pipe 
% input 
% L: pipe length (m),  R: pipe diameter (m),  K: parameter (N*s^n/m^2) 
% n: flow index,  delP: pressure drop (Pa) 
% output 
% avgv: avg. velocity (m/s) 
h = 2*R/100; r = [-R:h:R]; % h: step size 
Rn = R.^((n+1)/n); Pn = (delP./(2*K*L)).^(1/n); 
v = (1 - (abs(r)/R).^((n+1)/n)) .* Pn .* Rn * n/(n+1); 
avgv = Rn .* Pn .* n/(3*n+1); dvr = diff(v)/h; dvr = [dvr dvr(end)]; taurx = -K.*abs(dvr).^(n-1).*dvr;  
subplot(1,2,1), plot(v,r), grid, xlabel('v_x(r)'), ylabel('r'), axis tight 
subplot(1,2,2), plot(taurx,r), grid, xlabel('\tau_{rx}(r)'), ylabel('r'), 
axis tight 
end 
% slabtemp.m
LA = 0.015; LB = 0.1; LC = 0.075; Lt = LA+LB+LC; kA = 0.151; kC = 0.762;
qx = -15; T0 = 255; xspan = [0 Lt];
[x T] = ode45(@slabmT, xspan, T0, [], LA,LB,kA,kC,qx);
plot(x,T), grid, axis([0 Lt 250 310]), xlabel('x(m)'), ylabel('T(K)')

function dT = slabmT(x,T,LA,LB,kA,kC,qx) 
% slabmT.m: heat transfer through multilayer slab 
if x <= LA, dT = -qx/kA; 
elseif x <= LA+LB, dT = -qx/(2.5*exp(-1225/T)); 
else, dT = -qx/kC; end  
end 
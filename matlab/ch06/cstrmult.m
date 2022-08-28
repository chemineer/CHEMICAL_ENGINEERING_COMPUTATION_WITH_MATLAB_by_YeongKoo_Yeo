function fC = cstrmult(C,ka,kc,Ca0,Cb0,v0,V) 
% C(1)=Ca, C(2)=Cb, C(3)=Cb, C(4)=Cd 
ra = -ka*C(1)*C(2)^2 - 2*kc*C(1)^2*C(3)^3/3; rb = -2*ka*C(1)*C(2)^2; 
rc = ka*C(1)*C(2)^2 - kc*C(1)^2*C(3)^3; rd = kc*C(1)^2*C(3)^3/3; 
fC = [v0*Ca0 - v0*C(1) + ra*V; v0*Cb0 - v0*C(2) + rb*V; -v0*C(3) + rc*V; -v0*C(4) + rd*V]; 
end 
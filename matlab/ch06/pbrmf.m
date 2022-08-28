function dxdw = pbrmf(w,x,k,fa0,ca0,ka,kb,kc) 
kfc = k*ca0^2/fa0; 
dxdw = [kfc*(1-x(1))^2/(1+ka*ca0*(1-x(1)));          
kfc*(1-x(2))^2/(1+ka*ca0*(1-x(2))+kc*ca0*x(2));          
kfc*(1-x(3))^2/(1+ka*ca0*(1-x(3))+kb*ca0*(1-x(3)))^2;          
kfc*(1-x(4))^2/(1+ka*ca0*(1-x(4))+kb*ca0*(1-x(4))+kc*ca0*x(4))^2]; 
end 
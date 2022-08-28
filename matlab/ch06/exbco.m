function dz = exbco(V,z,Ca0,Fa0,Cp0,Cpc,Ua,m,T1,T2,k1,K2,E,R,dH,hx) 
% z(1)=Ta, z(2)=X, z(3)=T 
k = k1*exp(E*(1/T1 - 1/z(3))/R); Kc = K2*exp(dH*(1/T2 - 1/z(3))/R); ra = -k*Ca0*(1 - (1 + 1/Kc)*z(2)); 
if hx == 'co', dz(1) = Ua*(z(3) - z(1))/(m*Cpc); 
elseif hx == 'cn', dz(1) = -Ua*(z(3) - z(1))/(m*Cpc);  
else, dz(1) = 0; end 
dz(2) = -ra/Fa0; dz(3) = (ra*dH - Ua*(z(3)-z(1)))/(Fa0*Cp0); dz = dz'; 
end 
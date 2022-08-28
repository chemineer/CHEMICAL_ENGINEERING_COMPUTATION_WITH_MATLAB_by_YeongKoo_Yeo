function f = binflash(x,v,T,z,A,B,C) 
% v- fraction of feed vaporized, t: temp(deg.C), z: feed composition 
% A, B, C: parameters of Antoine equation  
% x(1): x1, x(2): y1, x(3): P 
for k = 1:2, Ps(k) = exp(A(k) - B(k)/(T + C(k))); end 
f = [x(1)*(1-v) + x(2)*v - z; x(1)*Ps(1) - x(2)*x(3); (1-x(1))*Ps(2) - (1-x(2))*x(3)]; 
end 
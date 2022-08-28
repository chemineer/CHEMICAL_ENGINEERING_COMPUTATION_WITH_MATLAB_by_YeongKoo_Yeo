function By = twophreg(Bx) 
% 2-phase flow regimes  
% Bx, By: Baker parameters 
By = []; 
C1 = exp(9.774459 - 0.6548*log(Bx)); 
C2 = exp(8.67964 - 0.1901*log(Bx)); 
C3 = exp(11.3976 - 0.6084*log(Bx) + 0.0779*log(Bx).^2); 
C4 = exp(10.7448 - 1.6265*log(Bx) + 0.2839*log(Bx).^2); 
C5 = exp(14.569802 - 1.0173*log(Bx)); 
C6 = exp(7.8206 - 0.2189*log(Bx)); 
By = [C1' C2' C3' C4' C5' C6']; 
end 
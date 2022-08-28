% Heat of Vaporization of Water 
T = 0:350; hv = hvapn(T,'H2O');  
plot(T,hv), xlabel('T(C)'), ylabel('Heat of vaporization of water(cal/g)'), grid on 
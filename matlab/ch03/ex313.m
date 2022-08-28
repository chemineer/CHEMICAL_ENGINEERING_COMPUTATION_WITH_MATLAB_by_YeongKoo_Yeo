% Thermal Conductivity of Water
T = [0:350]+273.15; k = condL(T,'water'); plot(T-273.15, k), xlabel('T(C)') 
ylabel('Thermal conductivity of water(\mucal/s/cm/C)'), grid on 
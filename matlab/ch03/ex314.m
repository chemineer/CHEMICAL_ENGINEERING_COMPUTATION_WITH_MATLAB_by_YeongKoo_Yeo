% Thermal Conductivity of Propane
T = [25:900]+273.15; k = condG(T,'propane'); plot(T-273.15, k), xlabel('T(C)'), axis tight, grid on 
ylabel('Thermal conductivity of propane(\mucal/s/cm/K)') 
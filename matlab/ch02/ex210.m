%Catalytic Dehydrogenation Reactions of Ethane C2H6 
K1  =  3.78; K2  =  0.137; x0  =  [0.80 0.10]; 
f  =  @(x) [x(1)*(x(1) + 2*x(2))/(1 - x(1) - x(2))/(1 + x(1) + 2*x(2)) - K1; x(2)*(x(1) + 2*x(2))^2/(1 - x(1) -x(2))/(1 + x(1) + 2*x(2))^2 - K2]; 
z  =  newtrapmv(f,x0)
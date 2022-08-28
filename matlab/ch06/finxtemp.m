% finxtemp.m
P = 162; Vspan = [0 4]; FN2 = 28.3; X0 = [38.3-FN2 0 0 1150]; pf = [P FN2];
[V X] = ode45(@adfun,Vspan,X0,[],pf);
xc = (X0(1) - X(:,1))/X0(1);
subplot(1,2,1), plot(V,X(:,4)), xlabel('Reactor volume(m^3)'), ylabel('Temperature(K)'), grid
subplot(1,2,2), plot(V,xc), xlabel('Reactor volume(m^3)'), ylabel('Conversion'), grid

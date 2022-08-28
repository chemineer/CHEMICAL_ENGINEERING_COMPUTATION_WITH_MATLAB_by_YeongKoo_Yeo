% convtemp.m
P = 162; Vspan = [0 4]; FN2 = [28.3 18.3 8.3 3.3 0.0]; nF = length(FN2);
for i = 1:nF
X0 = [38.3-FN2(i) 0 0 1150]; pf = [P FN2(i)];
[V X] = ode45(@adfun,Vspan,X0,[],pf);
xc(i) = (X0(1) - X(end,1))/X0(1); Tr(i) = X(end,4);
end
xc, Tr

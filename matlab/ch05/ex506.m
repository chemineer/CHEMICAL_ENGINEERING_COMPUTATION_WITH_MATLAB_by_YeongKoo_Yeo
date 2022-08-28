% ffcor.m: friction factor correlation 
Nre = [8500 20000 30000 60000 700000 1000000 10000000]; % Reynolds number 
f = [0.008 0.0065 0.006 0.005 0.003 0.0028 0.002]; % friction factor 
x = log10(Nre); y = log10(f); p = polyfit(x,y,1); % linear regression by polyfit 
a = 10^p(2), b = p(1) % identified parameters  
Nrex = linspace(min(Nre),max(Nre),100); fvx = a*Nrex.^b; % compute f by the correlation 
plot(Nre,f,'o',Nrex,fvx), legend('Data','Correlation') 
xlabel('N_{Re}(Reynolds number)'), ylabel('f(Friction factor)') 
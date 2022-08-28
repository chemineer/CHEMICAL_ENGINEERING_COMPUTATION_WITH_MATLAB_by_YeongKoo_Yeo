% rxntemfit.m: polynomial fitting by polyfit function 
t  =  [1 2 3 4 5 6 7 8]; tp  =  1:0.1:8; 
T  =  [50.8 56.4 55.1 60.6 61.5 59.5 54.1 53.8]; % temperature data 
p1  =  polyfit(t,T,1), p2  =  polyfit(t,T,2), p4  =  polyfit(t,T,4) % polynomial fitting 
T1  =  polyval(p1,tp); % Calculate temperature by 1st-order polynomial 
T2  =  polyval(p2,tp); % Calculate temperature by 2nd-order polynomial 
T4  =  polyval(p4,tp); % Calculate temperature by 4th-order polynomial 
plot(t,T,'o',tp,T1,':',tp,T2,'.-',tp,T4,'-'), xlabel('t(hr)'), ylabel('T(C)'), grid 
legend('Data','1st-order','2nd-order','4th-order')  
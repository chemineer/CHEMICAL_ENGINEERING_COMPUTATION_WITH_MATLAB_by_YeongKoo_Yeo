% fitrxndata.m: fit polynomials to the reaction data 
t  =  [0,2,5,6,13,20,24,30,35,41,50]; 
c  =  [0.86,0.61,0.47,0.39,0.25,0.18,0.15,0.12,0.10,0.09,0.08]; 
tmin  =  min(t); tmax  =  max(t); tint  =  linspace(tmin,tmax,100); 
g  =  {'--',':','.-'}; m  =  [3 4 5]; % 2nd-, 3rd-, and 4th-order polynomials 
for k  =  1:3, p  =  lspolfit(t,c,m(k)); cint{k}  =  polyval(p,tint); end 
% Plot data and fit polynomials 
plot(t,c,'o',tint,cint{1},g{1},tint,cint{2},g{2},tint,cint{3},g{3}) 
xlabel('t(min)'), ylabel('C_A(mol/liter)') 
legend('Data','2nd-order','3rd-order','4th-order','location','best') 
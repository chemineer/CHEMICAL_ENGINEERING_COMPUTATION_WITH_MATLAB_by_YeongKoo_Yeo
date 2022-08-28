% rootloc.m
Kc = [0:1:40];
for i = 1:length(Kc)
pol = [1 6 11 6+2*Kc(i)]; sol = roots(pol);
for j = 1:length(sol), r(i,j) = sol(j); end
end
plot(r,'*'), grid, ylabel('Imaginary part'), xlabel('Real part'), title('Root locus (Kc=1~40)')

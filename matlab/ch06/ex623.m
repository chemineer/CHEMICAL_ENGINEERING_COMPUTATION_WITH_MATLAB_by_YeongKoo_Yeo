% rxnLH.m: diffusion and reaction in a catalyst pellet by Langmuir-Hinshelwood kinetics 
clear all; 
n = 50; R = 0.2; k = 100; De = 0.25; Kr = 1e9; % data 
m = n+1; h = R/n; y0 = zeros(n+1,1); Ca0 = 5e-5*[1 2 3 4]; p = [];  
for i = 1:m, x(i) = h*(i-1); end % radial distance 
for i = 1:length(Ca0)    
y0 = Ca0(i)*ones(n+1,1); y = fsolve(@catLH,y0,[],De,n,h,k,Kr,Ca0(i)); p = [p y/Ca0(i)]; 
end 
x = x/R; plot(x,p(:,1),x,p(:,2),':',x,p(:,3),'--',x,p(:,4),'.-') 
xlabel('r/R'), ylabel('C_A/C_{A0}') 
legend('C_{A0} = 5x10^{-5} mol/cm^3','C_{A0} = 10x10^{-5} mol/cm^3',...    
'C_{A0} = 15x10^{-5} mol/cm^3','C_{A0} = 20x10^{-5} mol/cm^3','location','best') 
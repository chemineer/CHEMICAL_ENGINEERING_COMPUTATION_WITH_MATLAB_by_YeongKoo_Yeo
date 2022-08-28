% comfacng.m 
T = 60; P = 100:10:5000; nz = []; 
for k = 1:4, Sg = 0.5+(k-1)*0.1; nzv = ngasZ(T,P,Sg); nz = [nz nzv']; end 
plot(P,nz(:,1),P,nz(:,2),':',P,nz(:,3),'.-',P,nz(:,4),'--'), grid 
axis([100 5000 0 1.1]), legend('Sg=0.5','Sg=0.6','Sg=0.7','Sg=0.8','location','best')  
xlabel('P(psia)'), ylabel('Compressibility factor, Z')
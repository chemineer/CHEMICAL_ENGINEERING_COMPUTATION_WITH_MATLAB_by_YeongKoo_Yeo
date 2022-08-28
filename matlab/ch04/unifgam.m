function gam = unifgam(k,R,Q,nu,amn,Nc,x,T) 
% Estimation of activity coefficients using the UNIFAC method 
% input:  
% R,Q: vectors of volumes and surface areas for each functional group 
% nu: number of functional groups contained in each component 
% (row: functional groups k, column: components i) 
% amn: matrix of group interaction parameters 
% Nc: number of components 
% x: vector of liquid-phase mole fraction 
% T: temperature (K) 
% output: 
% gam: vector of activity coefficients for each component 
r = zeros(1,Nc); q = zeros(1,Nc); tau = zeros(k,k); ek = zeros(k,Nc);  
theta = zeros(1,k); beta = zeros(Nc,k); s = zeros(1,k); z = 10;  
for j = 1:Nc, r(j) = sum(R.*nu(:,j)'); q(j) = sum(Q.*nu(:,j)');  
end % row vector(r,q) 
for i = 1:k, ek(i,:) = nu(i,:)*Q(i)./q; end 
tau = exp(-amn/T);  
for i = 1:Nc    
for j = 1:k, beta(i,j) = sum(ek(:,i).*tau(:,j)); end 
end 
for i = 1:k, theta(i) = sum(x.*q.*ek(i,:))/sum(x.*q); end 
for i = 1:k, s(i) = sum(theta'.*tau(:,i)); end 
J = r/sum(r.*x); L = q/sum(q.*x); gamc = 1 - J + log(J) - 5*q.*(1 - J./L + log(J./ L)); gamr = []; 
for i = 1:Nc    
sumb = 0;     
for j = 1:k        
gamval = theta(j)*beta(i,j)/s(j) - ek(j,i)*log(beta(i,j)/s(j));         
sumb = sumb + gamval;    
end     
gamr = [gamr q(i)*(1 - sumb)]; 
end 
gam = exp(gamc + gamr); 
end 
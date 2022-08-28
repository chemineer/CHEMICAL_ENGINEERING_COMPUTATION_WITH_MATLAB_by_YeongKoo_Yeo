function frfactor(eD, Nre) 
% Friction factor correlations 
% inputs: 
% eD: surface roughness parameter 
% Nre: Reynolds number 
f0 = 5e-3; % initial guess 
fSC = 1./(log10(eD/3.7 - (5.02./Nre).*log10(eD/3.7+14.5/Nre))).^2 /16; 
% Shacham eqn. 

funCB = @(f) (1/sqrt(f) + 1.7372*log(eD/3.7 + 1.255/Nre/sqrt(f))); 
% Colebrook eqn. 
fCB = fzero(funCB, f0); 
funCBW = @(f) (1/sqrt(f) + 4*log10(eD + 4.67/Nre/sqrt(f)) - 2.28); 
% Colebrook- White eqn. 
fCBW = fzero(funCBW, f0); 
fHA = 1./(-3.6*log10(6.9/Nre + (eD/3.7)^(10/9)))^2; % Haaland eqn. 
Av = eD/3.7 + (6.7/Nre)^0.9; % Chen eqn. 
fCH = 1./(-4*log10(eD/3.7 - 5.02*log10(Av)/Nre)).^2; 
funNK = @(f) (1./sqrt(f) - 4*log10(Nre*sqrt(f)) + 0.4); % Nikuradse eqn. 
fNK = fzero(funNK, f0); 
fBS = 0.0791*Nre^(-1/4); % Blasius eqn. 
fprintf('Shacham eqn.: f = %g\n', fSC);  
fprintf('Colebrook eqn.: f = %g\n', fCB); 
fprintf('Colebrook-White eqn.: f = %g\n', fCBW); 
fprintf('Haaland eqn.: f = %g\n', fHA);  
fprintf('Chen eqn.: f = %g\n', fCH); 
fprintf('Nikuradse eqn.: f = %g\n', fNK);  
fprintf('Blasius eqn.: f = %g\n', fBS); 
end 
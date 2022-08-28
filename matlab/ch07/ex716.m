% SO2abs.m: absorption of SO2 by water
Ht = 0.6; H = 26; Gm = 206; Lm = 12240; x2 = 0.0; y1 = 0.03; y2 = 0.003; % data
Sf = Lm/(H*Gm); % absorption factor
Nt = Sf/(Sf-1)*log((1 - 1/Sf)*((y1 - H*x2)/(y2 - H*x2)) + 1/Sf); % number of transfer units
Hpd = Nt*Ht; % total packing height
fprintf('Number of theoretical transfer units = %g\n', Nt);
fprintf('Total packing height = %g m\n', Hpd);

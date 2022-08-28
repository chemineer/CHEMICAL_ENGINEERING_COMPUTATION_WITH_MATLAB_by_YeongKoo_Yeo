% datfis.m: creation of a FIS object with options 
indat = [rand(10,1) 10*rand(10,1)-5]; % input data (10x2) 
outdat = rand(10,1); % output data (10x1) 
opt = genfisOptions('GridPartition'); % create a FIS object using grid partitioning 
opt.NumMembershipFunctions = [3 5]; % 3 gaussian and 5 triangular mf’s for the 1st and 2nd input variable 
opt.InputMembershipFunctionType = ["gaussmf" "trimf"];  
exfis = genfis(indat,outdat,opt); % generate a FIS object 
[x1,mf1] = plotmf(exfis,'input',1); subplot(2,1,1), plot(x1,mf1), xlabel('input 1 (gaussmf)') 
[x2,mf2] = plotmf(exfis,'input',2); subplot(2,1,2), plot(x2,mf2), xlabel('input 2 (trimf)') 
% ethanrxn.m: ethane-steam cracking reaction 
x0  =  [0.001 0.001 0.001 0.993 1 0.0001 5.992 1 0.001 10 10 10]'; T  =  1000; 
%[x, fval] =  fsolve(@rxnfun, x0,[],T); % use fsolve to solve the equation system 
[x, fval] =  fsolve(@rxnfun, x0); % use fsolve to solve the equation system 
comp  =  {'CH4','C2H4','C2H2','CO2','CO','O2','H2','H2O','C2H6'}; 
lamda  =  {'lambda1','lambda2','lambda3'}; 
fprintf('\n i\tComp. \t Initial Val.\t\tFinal val.\n'); 
for k = 1:length(x)-3, fprintf('%d\t%s \t%12.9f\t%12.9f\n',k,comp{k},x0(k),x(k)); end 
fprintf('\n i\tLambda\tInitial Val.\tFinal val.\n'); 
for i =  k+1:length(x), fprintf('%g\t%s\t\t%4.1f\t%15.9f\n',i,lamda{i-length (comp)},x0(i),x(i)); end 
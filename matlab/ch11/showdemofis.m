% showdemofis.m
load demofis;
demofis.Rules % display all rules
figure(1), plotfis(demofis) % construct a block diagram representing the system figure(2)
subplot(2,2,1), plotmf(demofis,'input',1) % the membership function for the 1st input
subplot(2,2,2), plotmf(demofis,'input',2) % the membership function for the 2nd input
subplot(2,2,3), plotmf(demofis,'output',1) % the membership function for the 1st output
subplot(2,2,4), gensurf(demofis) % FIS outputs for input variables

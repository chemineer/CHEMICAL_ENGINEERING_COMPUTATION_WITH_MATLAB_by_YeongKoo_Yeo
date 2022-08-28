% mfplots.m 
x = 0:0.1:10;  
y1 = trimf(x,[3 6 8]); % Triangular  
subplot(2,2,1), plot(x,y1), xlabel('trimf, P = [3 6 8]'), ylim([-0.05 1.05]), title('Triangular') 
y2 = sigmf(x,[2 4]); % Sigmoid  
subplot(2,2,2), plot(x,y2), xlabel('sigmf, P = [2 4]'), ylim([-0.05 1.05]), title('Sigmoidal') 
y3 = gaussmf(x,[2 5]); % Gaussian  
subplot(2,2,3), plot(x,y3), xlabel('gaussmf, P = [2 5]'), title('Gaussian') 
y4 = smf(x,[1 8]); % S-membership 
subplot(2,2,4), plot(x,y4), xlabel('smf, P = [1 8]'), ylim([-0.05 1.05]), title('S- membership') 
% myffnet.m 
x = [4.99 5.15 5.30 5.45 5.59 5.73 5.89 6.07 6.24 6.39 6.53 6.65 6.74 6.82 6.90 6.97 7.03 7.09 7.15 7.20]; % input 
t = [7.31 7.43 7.57 7.71 7.85 7.99 8.11 8.18 8.14 8.01 7.76 7.46 7.16 6.84 6.52 6.20 5.89 5.57 5.26 4.95]; % target 
ffnet = feedforwardnet(10); % one hidden layer with 10 neurons 
ffnet = train(ffnet,x,t); % train the feedforward network 
view(ffnet) 
y = ffnet(x); 
perf = perform(ffnet,y,t)  % assess performance 
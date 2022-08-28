% testffnet.m 
clear all;   
x = [0 1.99 0.10 2.28 0.21 2.56 0.32 2.85 0.44 3.13 0.58 3.43 0.74 3.77...    
0.94 4.12 1.21 4.47 1.55 4.83]; % train data: input 
t = [5.05 9.47 5.66 8.97 6.32 8.46 6.94 8.01 7.57 7.65 8.21 7.38 8.84...    
7.19 9.43 7.10 9.90 7.10 9.98 7.22]; % train data: target
xv = [1.85 0.05 2.13 0.16 2.42 0.27 2.70 0.38 2.99 0.51 3.28 0.66 3.59...    
0.83 3.95 1.07 4.30 1.37 4.65 1.70]; % test data: input 
tv = [9.69 5.36 9.23 5.99 8.71 6.63 8.22 7.26 7.82 7.90 7.51 8.52 7.28...    
9.14 7.13 9.70 7.09 9.99 7.14 9.86]; % test data: target 
samnet = feedforwardnet(3,'trainlm'); % create a network having one hidden layer with 3 nodes 
samnet = train(samnet,x,t); % train the network 
ty = samnet(xv); % output from the trained network 
plot(xv,tv,'o',xv,ty,'*'), xlabel('Input'), ylabel('Target'), legend('Target data','Network output') 
perf = perform(samnet,ty,tv) % assess performance 
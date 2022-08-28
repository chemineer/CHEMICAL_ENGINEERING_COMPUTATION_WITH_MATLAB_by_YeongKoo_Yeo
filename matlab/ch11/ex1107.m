% regfit.m 
clear all;   
x = [1.85 0.05 2.13 0.16 2.42 2.56 0.32 0.38 2.99 0.51 3.28 0.66 0.74 3.77...    
0.94 4.12 4.30 1.37 4.65 1.70 0 1.99 0.10 2.28 0.21 0.27 2.70 2.85 0.44...    
3.13 0.58 3.43 3.59 0.83 3.95 1.07 1.21 4.47 1.55 4.83]; % train data: input 
t = [9.69 5.36 9.23 5.99 8.71 8.46 6.94 7.26 7.82 7.90 7.51 8.52 8.84 7.19...    
9.43 7.10 7.09 9.99 7.14 9.86 5.05 9.47 5.66 8.97 6.32 6.63 8.22 8.01...    
7.57 7.65 8.21 7.38 7.28 9.14 7.13 9.70 9.90 7.10 9.98 7.22]; % train data: target 
rgnet = feedforwardnet(10,'trainlm'); % create a network having one hidden layer with 10 nodes 
rgnet = train(rgnet,x,t); % train the network 
y = rgnet(x); % output from the trained network 
[r,m,b] = regression(t,y) % calculate regression between targets and outputs 
figure(1), plotregression(t,y) % plot the regression 
figure(2), plotfit(rgnet,x,t) % plot the output function
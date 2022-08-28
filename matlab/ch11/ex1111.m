% linsvmr.m 
n = 100; m = 50;  % dimension of input data set X 
nz = 0.1; rng(1); X = sprandn(n,m,nz); Y = 1.2*X(:,20) + 0.5*sin(X(:,40)) + 0.3*randn(n,1); 
Tmd = fitrlinear(X,Y, 'Holdout',0.3); % create SVM regression model (30% of data are used in testing) 
modl = Tmd.Trained{1}; 
indtrain = training(Tmd.Partition); indtest = test(Tmd.Partition); 
ytrain = predict(modl,X(indtrain,:)); ytest = predict(modl,X(indtest,:)); 
y = Y(indtest); % data 
k = 1: length(y); plot(k,ytest,k,y,'o'), xlabel('Sample no.'),ylabel('Y'),legend('Estimation','Data') 
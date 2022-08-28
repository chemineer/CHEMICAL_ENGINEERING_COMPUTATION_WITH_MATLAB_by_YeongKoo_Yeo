% estmoist.m: predict output by SVM regression model 
clear all; D = datgen; % generate data set 
[n,m] = size(D); % data set size (n: number of rows, m: number of columns)  
X = [D(:,1) D(:,2) D(:,3)]; Y = D(:,4); 
rng(10); % specify positive integer seed value 
testdat = cvpartition(n,'Holdout',0.3); % hold out 30% of the data for testing 
indtrain = training(testdat); % training set indices 
indtest = test(testdat); % test set indices 
modl = fitrsvm(X(indtrain,:),Y(indtrain),'Standardize',true); % generation of regression model
yp = predict(modl,X(indtest,:)); % predict outputs  
yd = D(indtest,4); k = 1:length(yd);  
plot(k,yp,k,yd,'o')  % plot data and estimation results 
xlabel('Sample no.'),ylabel('Moisture(%)'),legend('Estimation','Data','Location','best') 
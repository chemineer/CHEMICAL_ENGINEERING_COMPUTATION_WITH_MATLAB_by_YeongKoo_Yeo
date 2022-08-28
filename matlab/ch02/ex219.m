% ttestex.m: perform t-test 
abdat  =  [0.72, 0.54, 0.62, 0.80, 0.76, 0.64, 0.75, 0.94, 0.85, 0.44]; 
avg  =  mean(abdat); stdv  =  std(abdat); % mean and standard deviation 
mu  =  0.86; % null hypothesis 
cf  =  0.05; % significance level 
%[h,p]  =  ttest(abdat,mu,cf,'both'); % both: both tails
[h,p]  =  ttest(abdat,mu,'tail','both'); % both: both tails 
fprintf('h (result of the hypothesis test)  =  %g\n',h); 
fprintf('p (probability)  =  %g\n',p);  
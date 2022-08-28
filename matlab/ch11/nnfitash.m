% nnfitash.m 
clear all;  
xc1 = [0.0015 0.0076 0.0214 0.0273 0.0304 0.0263 0.0131 0.0087 0.0097 0.0069...    
0.0015 0.0143 0.0133 0.0083 0.0047 0.0039 0.0043 0.0032 0.0006 0.0056...    
0.0169 0.0240 0.0208 0.0125 0.0077 0.0086 0.0097 0.0117 0.0127 0.0083...    
0.0123 0.0243 0.0555 0.0619 0.0666 0.0763 0.0796 0.0820 0.0901 0.1051...    
0.1232 0.1313 0.1319 0.1387 0.1531 0.1636 0.1766 0.1940 0.2035 0.2092]; 
xc2 = [0.0365 0.0363 0.0361 0.0359 0.0357 0.0355 0.0353 0.0351 0.0350 0.0319...    
0.0169 0.0121 0.0119 0.0160 0.0227 0.0226 0.0224 0.0222 0.0220 0.0218...    
0.0216 0.0214 0.0192 0.0137 0.0149 0.0225 0.0224 0.0195 0.0120 0.0045...    
0.0009 0.0013 0.0024 0.0027 0.0030 0.0033 0.0036 0.0038 0.0041 0.0044...    
0.0046 0.0049 0.0052 0.0055 0.0057 0.0060 0.0063 0.0065 0.0068 0.0071]; 
xc3 = [0.0003 0.0007 0.0010 0.0013 0.0017 0.0020 0.0024 0.0027 0.0031 0.0034...    
0.0041 0.0044 0.0048 0.0051 0.0054 0.0058 0.0061 0.0065 0.0068 0.0072...    
0.0075 0.0078 0.0082 0.0085 0.0089 0.0092 0.0095 0.0099 0.0102 0.0106...    
0.0109 0.0112 0.0119 0.0123 0.0126 0.0129 0.0133 0.0136 0.0140 0.0143... 
0.0147 0.0150 0.0153 0.0157 0.0161 0.0160 0.0167 0.0170 0.0174 0.0177]; 
t = [0.1113 0.1131 0.1173 0.1182 0.1169 0.1190 0.1227 0.1257 0.1327 0.1350...    
0.1346 0.1327 0.1308 0.1302 0.1310 0.1324 0.1330 0.1318 0.1316 0.1320...    
0.1305 0.1306 0.1277 0.1254 0.1256 0.1230 0.1225 0.1207 0.1223 0.1279... 
0.1316 0.1366 0.1373 0.1379 0.1369 0.1374 0.1331 0.1303 0.1239 0.1174...    
0.1135 0.1120 0.1095 0.1078 0.1050 0.1016 0.0964 0.0940 0.0885 0.0871]; 
x = [xc1; xc2; xc3]; % input data matrix 
ppnet = fitnet(10); % one hidden layer with 10 neurons 
ppnet.divideParam.trainRatio = 0.7; % use 70% for training 
ppnet.divideParam.valRatio = 0.15; % use 15% for validation 
ppnet.divideParam.testRatio = 0.15; % use 15% for testing 
[ppnet,tr] = train(ppnet,x,t); 
y = ppnet(x); % output calculation by network 
neterr = gsubtract(y,t); % error 
perm = perform(ppnet,t,y) % performance 
view(ppnet) 
figure(1), plotperform(tr) 
figure(2), plotfit(ppnet,t,y) 
figure(3), plotregression(t,y) 
figure(4), ploterrhist(neterr) 
% Compare data and network output 
testx = x(:,tr.testInd); % input data set for test 
testt = t(:,tr.testInd); % target data set for test 
testy = ppnet(testx); % network output for test input data 
k = 1:length(tr.testInd);  
figure(5), plot(k,testy,k,testt,'o') 
xlabel('Sample no.'),ylabel('Ash(%)'),legend('Estimation by network','Data','Location','best')  
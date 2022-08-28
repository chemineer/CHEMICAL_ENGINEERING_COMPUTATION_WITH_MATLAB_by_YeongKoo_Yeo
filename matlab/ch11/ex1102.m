% fzopn.m: fuzzy arithmetic operations using fuzarith function
N = 101; minx = -20; maxx = 20;
x = linspace(minx,maxx,N); % generate N points between minx and maxx
A = trapmf(x,[-10 -2 1 3]); B = gaussmf(x,[2 5]); % define A and B using membership functions
Csum = fuzarith(x,A,B,'sum'); % A+B
Csub = fuzarith(x,A,B,'sub'); % A-B
Cprod = fuzarith(x,A,B,'prod'); % A*B
subplot(2,2,1), plot(x,A,'b--',x,B,'m:',x,Csum,'k-'), title(' A+B'), legend('A','B','A+B')
subplot(2,2,2), plot(x,A,'b--',x,B,'m:',x,Csub,'k-'), title(' A-B'), legend('A','B','A-B')
subplot(2,2,3), plot(x,A,'b--',x,B,'m:',x,Cprod,'k-'), title(' A*B'), legend('A','B','A*B')


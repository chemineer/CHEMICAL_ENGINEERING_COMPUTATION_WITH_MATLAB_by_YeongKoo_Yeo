% plotnyquist
w = logspace(-2, 1, 300); num = 2; den = conv([10 1], [2.5 1]); [x, p] = bode(num,den, w);
subplot(1,2,1), polar((pi/180)*p, x), title('polar plot of a 2nd-order process')
subplot(1,2,2), nyquist(num, den)

% highresp.m
num = 1; den1 = conv([2 1],[2 1]); den2 = conv(den1,den1); den3 = conv(den2,[2 1]); t = [0:0.1:20];
y1 = step(num, den1, t); y2 = step(num, den2, t); y3 = step(num, den3, t);
plot(t, y1, ':', t, y2, '--', t, y3), grid, xlabel('Time t(sec)'), ylabel('Response y(t)')
title('Step response of higher-order process (n=2, 4, 5)'), legend('n = 2', 'n = 4', 'n = 5', 'location','best')

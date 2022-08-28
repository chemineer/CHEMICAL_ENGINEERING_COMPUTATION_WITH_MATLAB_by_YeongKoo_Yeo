% levelcon.m
Kc1 = 5; Kc2 = 20; Kc3 = 50; t = [0:0.1:10]; num1 = Kc1; num2 = Kc2; num3 = Kc3;
den1 = [5 1+Kc1]; den2 = [5 1+Kc2]; den3 = [5 1+Kc3];
y1 = step(num1, den1, t); y2 = step(num2, den2, t); y3 = step(num3, den3, t);
plot(t, y1, ':', t, y2, '--', t, y3), xlabel('Time t(sec)'), ylabel('output y (t)'), legend('Kc = 5', 'Kc = 20', 'Kc = 50')

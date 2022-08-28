% simstep.m
num = 3; den = [2 1]; t = [0:0.1:10];
y = step(num, den, t);
plot(t, y), grid, title('Step response of a 1st-order process'), xlabel('t (sec)'), ylabel('Response y(t)')

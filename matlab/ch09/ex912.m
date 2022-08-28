% step2ndpro.m
tau = 0.5; z1 = 0.5; z2 = 1; z3 = 1.5; num = 1; t = [0:0.1:10];
den1 = [tau*tau 2*tau*z1 1]; den2 = [tau*tau 2*tau*z2 1]; den3 = [tau*tau 2*tau*z3 1];
y1 = step(num, den1, t); y2 = step(num, den2, t); y3 = step(num, den3, t);
plot(t, y1, ':', t, y2, '--', t, y3), grid, xlabel('Time t(sec)'), ylabel('Output y(t)')
title('Step responses of a 2nd-order process');
legend('Underdamped', 'Critically damped', 'Overdamped');

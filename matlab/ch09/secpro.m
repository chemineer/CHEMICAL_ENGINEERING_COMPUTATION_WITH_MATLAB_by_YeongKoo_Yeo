% secpro.m : feedback control of 2nd-order process using P-controller
% Input controller gain
Kc = input('Controller gain = '); num = 0.5*Kc; den = [0.5 1 0.5*Kc];
% Unit step response
t = [0:0.1:10]; y = step(num, den, t);
plot(t, y), grid, xlabel('Time t(sec)'), ylabel('Output C(t)');

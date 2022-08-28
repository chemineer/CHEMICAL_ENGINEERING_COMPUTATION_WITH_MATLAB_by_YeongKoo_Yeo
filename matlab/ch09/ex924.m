% delaybode.m
num = [0.4 1]; den = conv([0.3 1], conv([1 1], [1 1]));
G = tf(num, den, 'iodelay', 0.2); [mag, phase, w] = bode(G);
subplot(2,1,1), loglog(w, squeeze(mag)), grid, ylabel('Amplitude'), xlabel('Frequency(rad/time)')
subplot(2,1,2), semilogx(w, squeeze(phase)), grid, ylabel('Phase(deg)'),
xlabel('Frequency(rad/time)')

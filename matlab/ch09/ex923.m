% bodezeta.m
w = logspace(-2, 1, 300); zeta = [0:0.25:1]; num = 1;
for i = 1:length(zeta),
den = [2.25 3*zeta(i) 1]; [ar, phase] = bode(num, den, w);
lar = 20*log10(ar);
subplot(2,1,1), loglog(w,ar), xlabel('w(rad/min)'), ylabel('AR'), title('Response of 2nd-order process')
text(0.8, 20, 'zeta=0'), text(0.4, 0.2,'zeta=1'), grid, hold on
subplot(2,1,2), semilogx(w, phase), grid, xlabel('w(rad/min)'), ylabel('\phi'), hold on
end
hold off

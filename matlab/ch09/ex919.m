% tunpid.m
tau = 2; zeta = 0.25; Kp = 3; Kc = 5; h(1,:) = '- '; h(2,:) = ': '; h(3,:) = '-.'; h(4,:) = '--';
t0 = 0; delt = 0.1; fint = 20; ms = 1;
% Constant reset time (tauI)
tauI = 1; tauD = [0.5 1 5 10]; % try 4 different tauD

subplot(1,2,1)
for i = 1:length(tauD)
num = Kc*Kp*[tauI*tauD(i) tauI 1]; d1 = tauI*tau^2;
d2 = 2*tauI*tau*zeta+Kc*Kp*tauI*tauD(i); d3 = tauI*(1+Kc*Kp); d4 = Kc*Kp;
den = [d1 d2 d3 d4];
[y,t] = stepnp(num, den, t0, delt, fint, ms); plot(t,y,h(i,:)), hold on
end
st = 1+0*t; plot(t,st), hold off, legend('\tau_D=0.5','\tau_D=1','\tau_D=5','\tau_D=10','location','best')
xlabel('Time(min)'), ylabel('Output, y(t)'), title('PID control(Kc=5,\tau_I=1)')
% Constant derivative time (tauD)
tauD = 0.5; tauI = [0.5 1 5 10]; % try 4 different tauI
subplot(1,2,2)
for i = 1:length(tauI)
num = Kc*Kp*[tauI(i)*tauD tauI(i) 1]; d1 = tauI(i)*tau^2;
d2 = 2*tauI(i)*tau*zeta+Kc*Kp*tauI(i)*tauD;
d3 = tauI(i)*(1+Kc*Kp); d4 = Kc*Kp; den = [d1 d2 d3 d4];
[y,t] = stepnp(num, den, t0, delt, fint, ms); plot(t,y,h(i,:)), hold on
end
st = 1+0*t; plot(t,st), hold off, legend('\tau_I=0.5','\tau_I=1','\tau_I=5','\tau_I=10','location','best')
xlabel('Time(min)'), ylabel('Output, y(t)'), title('PID control(Kc=5,\tau_D=0.5)')

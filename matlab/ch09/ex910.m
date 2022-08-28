% sin1st.m: sinusoidal response of the 1st-order process
num = 3; den = [2 1]; G = tf(num,den);
t = [0:0.1:10]; u = sin(3*t); z = t*0; y = lsim(G,u,t);
plot(t,y,t,u,':'), legend('Output y(t)','Input u(t)'), hold on
plot(t,z), hold off, xlabel('t(sec)'), ylabel('Response y(t)')
title('Sinusoidal response of 1st-order process')

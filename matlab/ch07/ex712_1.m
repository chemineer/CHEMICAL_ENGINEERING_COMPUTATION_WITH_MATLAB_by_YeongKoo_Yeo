% absbz.m
A = [-0.325 0.125 0 0 0;0.2 -0.325 0.125 0 0;0 0.2 -0.325 0.125 0; 0 0 0.2 -0.325 0.125;0 0 0 0.2 -0.325];
B = [0.2 0;0 0;0 0;0 0;0 0.25]; C = [0 0 0 0 1;0.5 0 0 0 0]; D = [0 0;0 0];
Us = [0.0;0.1]; xs = -inv(A)*B*Us; ys = C*xs + D*Us; % Steady-state
[y,x,t] = step(A,B,C,D,2); y = 0.05*y; x = 0.05*x;
for k = 1:2, y(:,k) = y(:,k) + ys(k); end; % Y(t) = ys + y(t)
for k = 1:5, x(:,k) = x(:,k) + xs(k); end; % X(t) = xs + x(t)
subplot(2,2,1), plot(t,y(:,1)), xlabel('t(min)'), ylabel('x_5'), grid, axis tight
subplot(2,2,2), plot(t,y(:,2)), xlabel('t(min)'), ylabel('y_1'), grid, axis tight
subplot(2,2,3), plot(t,x(:,1),t,x(:,2),'--',t,x(:,3),'.-',t,x(:,4),'*-',t,x(:,5),'o-')
xlabel('t(min)'), ylabel('x'), axis tight, legend('x_1','x_2','x_3','x_4', 'x_5')

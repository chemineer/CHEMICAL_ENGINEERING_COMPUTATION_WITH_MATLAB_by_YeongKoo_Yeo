% catwt.m 
ca0 = 1; fa0 = 1.5; k = 10; ka = 1; kb = 2; kc = 20; wf = 2;  % data 
wint = [0 wf]; x0 = [0 0 0 0]; 
[w x] = ode45(@pbrmf,wint,x0,[],k,fa0,ca0,ka,kb,kc); 
plot(w,x(:,1),w,x(:,2),':',w,x(:,3),'.-',w,x(:,4),'--'), xlabel('W(kg)'), ylabel('X') 
legend('X_1','X_2','X_3','X_4','location','best') 
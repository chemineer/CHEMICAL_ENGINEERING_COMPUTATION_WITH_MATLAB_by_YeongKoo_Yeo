% Use of trapz and cumtrapz
t  =  [0 0.5 1.2 1.6 2.5 3.1 4.8 6.9]; t1  =  [0 0.5 1.2 1.6 2.5];
v  =  [0 5 12 15 23 28 38 47]; v1  =  [0 5 12 15 23];
z1  =  trapz(t1,v1)
z  =  cumtrapz(t,v)
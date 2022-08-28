ng = [2 1]; dg = [1 3 2]; nh = 1; dh = [1 1]; [nt, dt] = feedback(ng, dg, nh, dh, -1);
Gcl = tf(nt, dt)

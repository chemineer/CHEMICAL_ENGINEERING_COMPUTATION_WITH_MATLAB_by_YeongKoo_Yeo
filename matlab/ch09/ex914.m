% tdelay.m
G1 = tf(3, [3 1]); delay = 1.6; set(G1, 'iodelay', delay);
G2 = tf(3, conv(conv([0.1 1], [0.5 1]), conv([1 1], [3 1])));
step(G1); hold on; step(G2, ':'); legend('G1', 'G2'), grid, hold off

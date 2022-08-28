G = tf(0.8, [5 1]); delay = 2; [m, p, w] = bode(G); Mag = m(1, :);
Phase = p(1, :) - ((180/pi)*delay*w'); [Gm, Pm, Wcg, Wcp] = margin(Mag, Phase, w)
G = tf(0.8, [5 1]); theta = 2; set(G, 'iodelay', theta); margin(G)

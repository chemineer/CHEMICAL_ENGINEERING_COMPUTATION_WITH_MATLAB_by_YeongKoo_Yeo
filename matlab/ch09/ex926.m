G = tf(12.76,[5 1]); set(G,'iodelay',1); nyquist(G), xlabel('Real axis'), ylabel('Imaginary axis')
title('Nyquist diagram of a time delay')

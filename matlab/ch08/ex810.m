% htdiri.m
alpha = 4.52e-7; nx = 25; dx = 0.5/nx; dt = 300; nt = 18000/dt;
u0 = 100*ones(1,nx+1); bci = 18; bcf = 18;
[u r] = parabDbc(nx, nt, dx, dt, alpha, u0, bci, bcf);
surf(u), axis([0 nx+1 0 nt+1 0 110]), view([-217 30]), xlabel('x'), ylabel('t'), zlabel('T'), r

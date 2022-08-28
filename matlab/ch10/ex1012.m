% optcrude.m: LP by interior point method
A = zeros(8,13);
A(1:3,1:5) = [0.6 0.5 0.3 0.4 0.4;0.2 0.2 0.3 0.3 0.1;0.1 0.2 0.3 0.2 0.2];
A(1:3,6:8) = eye(3); A(4:8,1:5) = eye(5); A(4:8,9:13) = eye(5);
b = 1e4*[17 8.5 7.5 8 10 10 10 6]'; c = [-31.1 -20.7 -21.3 -23.2 -20.2]; tol = 1e-6;
[xopt,fopt,basic] = barnslp(A,b,c,tol)

% c5distBPeu.m: calculation of 5-component distillation using BP method
% (English unit)
% Operating conditions and parameters
clear all;
opdat.N = 16; opdat.nc = 5; % N: total stages, nc: number of components
intv = zeros(1,opdat.N); intc = zeros(opdat.nc,opdat.N);
opdat.F = intv; opdat.Tf = intv; opdat.Pf = intv; opdat.hF = intv;
opdat.P = intv; opdat.T = intv; opdat.V = intv; opdat.L = intv;
opdat.U = intv; opdat.W = intv; opdat.Q = intv; opdat.eos = 'pr';
opdat.hV = intv; opdat.hL = intv; opdat.nf = 0; opdat.x = intc; opdat.y = intc;
opdat.z = intc;
% Physical properties and parameters for each component
mxdat.Pc = [709.8 617.4 550.7 489.5 440.0]; % critical pressure (psia)
mxdat.Tc = [550.0 665.9 765.3 845.9 914.2]; % critical temperature (R)
mxdat.k = zeros(opdat.nc,opdat.nc); % matrix of binary interaction parameters
mxdat.w = [0.1064 0.1538 0.1954 0.2387 0.2972]; % acentric factors
% Antoine equation parameters
mxdat.Ant = [5.38389 2847.921 434.898; 5.35342 3371.084 414.488;
5.74162 4126.385 409.5179; 5.853654 4598.287 394.4148; 6.03924 5085.758 382.794];
% Specific heat parameters (a(i))
mxdat.Afi = [11.51606 0.140309e-1 0.0854034e-4 -0.110608e-7 0.316220e-11;
15.58683 0.2504953e-1 0.1404258e-4 -0.352626e-7 1.864467e-11;
20.79783 0.314329e-1 0.192851e-4 -0.458865e-7 2.380972e-11;
25.64627 0.389176e-1 0.239729e-4 -0.584262e-7 3.079918e-11
30.17847 0.519926e-1 0.030488e-4 -0.27640e-7 1.346731e-11];
% Operating conditions (P: psia, T: F, flow: lbmol/hr)
opdat.nf = [6 9]; % feed stage
opdat.fstate = 'V'; % feed state
opdat.F(opdat.nf) = [41 59]; opdat.V(1) = 15; opdat.L(1) = 150;
opdat.U(1) = 5; opdat.U(3) = 3; opdat.W(13) = 37;
opdat.L(opdat.N) = sum(opdat.F(opdat.nf)) - opdat.V(1) - opdat.U(1) - opdat.U(3) - opdat.W(13);
opdat.Pf(opdat.nf) = [300 275]; opdat.Tf(opdat.nf) = [170 230];
opdat.P = 240*ones(1,opdat.N); opdat.Q(3) = 2e5;
opdat.z(:,opdat.nf(1)) = [0.061 0.342 0.463 0.122 0.012];
opdat.z(:,opdat.nf(2)) = [0.0085 0.1017 0.3051 0.5085 0.0762];
% Initialization and guess
opdat.T0 = (300/opdat.N)*[1:opdat.N]; % guess stage temperatures (F)
opdat.V(2:opdat.N-1) = 170; opdat.criv = 1e-3; % convergence criterion
% BP calculation
[x,y,T,L,V,iter] = distBPeu(opdat,mxdat);
fprintf('Number of iterations (convergence criterion: %g): %d\n', opdat.criv, iter)
disp('Temperature: '), T
Nx = 1:opdat.N;
x1 = x(1,:); x2 = x(2,:); x3 = x(3,:); x4 = x(4,:); x5 = x(5,:);
y1 = y(1,:); y2 = y(2,:); y3 = y(3,:); y4 = y(4,:); y5 = y(5,:);
subplot(2,2,1), plot(Nx,T), xlabel('Stage'),ylabel('T(F)'), axis tight
subplot(2,2,2), plot(x1,Nx,x2,Nx,'--',x3,Nx,'.-',x4,Nx,':',x5,Nx,'*-')
xlabel('x'),ylabel('Stage'), legend('x_1','x_2','x_3','x_4','x_5'), axis
([0 1 1 16])
subplot(2,2,3), plot(y1,Nx,y2,Nx,'--',y3,Nx,'.-',y4,Nx,':',y5,Nx,'*-')
xlabel('y'),ylabel('Stage'), legend('y_1','y_2','y_3','y_4','y_5'), axis
([0 1 1 16])

% pintemp.m: pin temperature profile
Ta = 25; L = 0.01; Ac = 6.28e-5; Ak = 3.14e-6; hc = 20; % data
rho = 2707; k = 220; Cp = 896; % properties of pure aluminum
V = Ak*L; Ch = rho*Cp*V; Rc = 1/(hc*Ac); Rk = L/(k*Ak);
T0 = Ta; Tb = 100; % initial pin and base temperature at t = 0
tspan = [0 8]; dT = @(t,T) (1/(Ch*Rc))*Ta + (1/(Ch*Rk))*Tb - (1/(Ch*Rc) + 1/ (Ch*Rk))*T;
[t,T] = ode45(dT,tspan,T0); plot(t,T), xlabel('t(s)'), ylabel('T(deg.C)')
Ts = ((1/(Ch*Rc))*Ta + (1/(Ch*Rk))*Tb)/(1/(Ch*Rc) + 1/(Ch*Rk)); % st-st temp.
fprintf('At steady-state, T(pin) = %g deg.C\n',Ts);

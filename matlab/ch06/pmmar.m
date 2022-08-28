function df = pmmar(t,x,pm)
% The reaction model is based on the following references:
% [1] Seth,V. and Gupta,S.K., J. of Polymer Eng, 15(3-4), pp. 283-323 (1995).
% [2] Ghosh,P., Gupta,S.K. and Saraf,D.N., Chemical Eng. J., 70, pp. 25-35 (1998).
% [3] Ray,A.B., Saraf,D.N. and Gupta,S.K., Polymer Eng. & Science, 35, pp. 1290- 1299 (1995).
% Retrieve data
M0 = pm.M0; MWm = pm.MWm; MWi = pm.MWi; Mjp = pm.Mjp; rhop = pm.rhop;
Vms = pm.Vms; Vps = pm.Vps; Vis = pm.Vis; T = pm.T; gam = 1; eff0 = 1;
% Assign variables
I = x(1); M = x(2); R = x(3); L0 = x(4); L1 = x(5); L2 = x(6);
N0 = x(7); N1 = x(8); N2 = x(9); Qg = x(10);
% Set values of temperature dependent parameters
rhom = 966.5 - 1.1*(T - 273.1); % monomer density (kg/m^3)
kd = 1.69e14 * exp(-125400/(8.314*T)); kp0 = 491.7 * exp(-18220/(8.314*T));
ktd0 = 9.8e4 * exp(-2937/(8.314*T));
Vm = 0.149 + 2.9e-4 * (T - 273.1); Vp = 0.0194 + 1.3e-4 * (T - 273.1 - 105);
% Ref. [3]
thet = 10^(124.1 - 1.0314e5 / T + 2.2735e7 / T^2); % theta_t
thep = 10^(80.3 - 7.5e4 / T + 1.765e7 / T^2); % theta_p
thef = 1e-3 * 10^(-40.86951 + 1.7179e4 / T); % theta_f
% Set parameters (Ref.[1])
gv = gam/Vp; eta13 = Vms*MWm/(Vps*Mjp); etai3 = Vis*MWi/(Vps*Mjp);
V = M*MWm/rhom + (M0 - M)*MWm/rhop; % volume of mixture (m^3)
psim = M*MWm/(rhom*V); psip = 1 - psim;
if L0+N0 == 0, rm = 0; else rm = (L1+N1)/(L0+N0); end % Ref. [1]
% Ref. [3]
ps = gam*(rhom*psim*Vms/eta13 + rhop*psip*Vps)/(rhom*psim*Vms*Vm + rhop*psip*Vps*Vp);
eff = eff0/(1 + thef*(M/V)/exp(etai3*(-ps+gv)));
% Reaction rate constants
ktd = 1/(1/ktd0 + thet*rm^2*(L0/V)/exp(-ps+gv));
kp = 1/(1/kp0 + thep*(L0/V)/exp(eta13*(-ps+gv)));
% Ref. [2]
ki = kp; kf = 0; ktc = 0;
if t<60, fr = M0*0.01*100/(242*60);  % Ref. [1]
else fr = 0; end
% Define differential equations (Ref. [1])
df(1,1) = -kd*I + fr;
df(2,1) = -(kp+kf)*L0*M/V - ki*R*M/V;
df(3,1) = 2*eff*kd*I - ki*R*M/V;
df(4,1) = ki*R*M/V - ktd*L0^2/V;
df(5,1) = ki*R*M/V + kp*M*L0/V - ktd*L0*L1/V + kf*M*(L0-L1)/V;
df(6,1) = ki*R*M/V + kp*M*(L0+2*L1)/V - ktd*L0*L2/V + kf*M*(L0-L2)/V;
df(7,1) = kf*M*L0/V + (ktd+ktc/2)*L0^2/V;
df(8,1) = kf*M*L1/V + ktd*L0*L1/V;
df(9,1) = kf*M*L2/V + ktd*L0*L2/V + ktc*L1^2/V;
df(10,1) = -57700*(-(kp+kf)*L0*M/V - ki*R*M/V);
end

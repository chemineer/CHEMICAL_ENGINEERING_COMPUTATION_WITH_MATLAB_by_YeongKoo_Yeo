%crfdata.m: cross-flow calculation data
t = 0.00254; % membrane thickness (cm)
Pm = [50 5]*1e-10; % permeability(cm^3*cm/(s*cm^2*cmHg)
alpa = Pm(1)/Pm(2); % ratio of permeabilities
ph = 80; % feed side pressure(cmHg)
pl = 20; % permeate side pressure(cmHg)
r = pl/ph; % pressure ratio (Plow/Phigh)
qf = 1e4; % feed rate(cm^3/s(STP))
xf = 0.5; % Feed composition (mole fraction)
theta = []; % stage-cut
xr = 0.25; % desired reject composition (mole fraction)

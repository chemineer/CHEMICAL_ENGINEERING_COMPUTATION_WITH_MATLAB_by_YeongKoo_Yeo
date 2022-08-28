% cpmdata.m: gas mixture and membrane data
t = 0.00254; % membrane thickness (cm)
Pm = [50 5]*1e-10; % permeability (cm^3*cm/(sec*cm^2*cmHg)
alpa = Pm(1)/Pm(2);
ph = 80; % feed-side pressure (cmHg)
pl = 20; % permeate-side pressure (cmHg)
r = pl/ph; % pressure ratio (Plow/Phigh)
qf = 1e4; % feed flow rate (cm^3/sec(STP))
xf = 0.5; % feed composition (mole fraction)
xr = 0.25; % desired reject composition (mole fraction)

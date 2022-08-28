% hxndat.m: Shell-and-tube heat exchanger data
% Two values of reference temperatures and physical properties at the two
% reference temperatures for tube and shell fluids
Trt = [323 283]; Trs = [375 289]; % reference temperatures (tube ans shell)(K)
rhoreft = [988.1 999.7]; rhorefs = [798 885]; % densities (kg/m^3)
mureft = [0.6 1.26]; murefs = [0.258 0.679]; % viscosities (mNs/m^2)
xkreft = [0.64 0.603]; xkrefs = [0.126 0.163]; % thermal conductivities
cpreft = [4183 4195]; cprefs = [1980 1675]; % heat capacities (J/kg/K)
% Heat exchanger geometry
Do = 25.4; Di = 19.86; % tube outside diameter (Do,mm) and inside diameter (Di,mm)
Xkw = 45; % heat conductivity of tube wall (W/m/K)
L = 2; % tube length (m)
Ls = 27; % tube sheet thickness (mm)
rf = 0.025; % tube roughness(mm)
Nt = 86; % total number of tubes in tube bundle
Layout = 1; % tube layout (1:triangular, 2:in-line square, 3:rotated square)
Pt = 1.25*Do; % tube pitch (mm)
Nss = 0; % number of pairs of sealing strips
Npass = 4; % number of passes
Rdt = 0.00036; Rds = 0.00018; % fouling resistance (m^2*K/W)
Ds = 305; % shell inside diameter(mm)
Dotl = 294; % shell outside tube limit(mm)
Dsb = 4.45; % shell-baffle clearance(mm)
Lbin = 165; Lbout = 165; % inlet and outlet baffle spacing(mm)
Lbc = 450; % central baffle spacing(mm)
Lc = 0.25*Ds; % baffle cut(mm)


% Specification of key variables
Ti1 = 298; Ti2 = 303; % inlet and outlet temperatures for tube-side(K)
Ts1 = 353; % Shell-side(hot stream) inlet temperature(K)
Ws = 4.8; % Shell-side(hot stream) mass flow rate(kg/sec)
fsT = 1; fsS = 1; % state of tube-side and shell-side fluids(1:liquid, 2:vapor)
ptype = 8; % problem type
Ts2 = 338; % guess outlet temperature of shell-side fluid
Wi = 5; % guess tube-side flow rate

class HeatExchangerData:
    def __init__(self):
        # Two values of reference temperatures and physical properties at the two
        # reference temperatures for tube and shell fluids
        self.Trt = [323, 283]  # reference temperatures (tube) (K)
        self.Trs = [375, 289]  # reference temperatures (shell) (K)
        
        self.rhoreft = [988.1, 999.7]  # densities (kg/m^3)
        self.rhorefs = [798, 885]      # densities (kg/m^3)
        
        self.mureft = [0.6, 1.26]      # viscosities (mNs/m^2)
        self.murefs = [0.258, 0.679]   # viscosities (mNs/m^2)
        
        self.xkreft = [0.64, 0.603]    # thermal conductivities
        self.xkrefs = [0.126, 0.163]   # thermal conductivities
        
        self.cpreft = [4183, 4195]     # heat capacities (J/kg/K)
        self.cprefs = [1980, 1675]     # heat capacities (J/kg/K)
        
        # Heat exchanger geometry
        self.Do = 25.4   # tube outside diameter (Do, mm)
        self.Di = 19.86  # tube inside diameter (Di, mm)
        self.Xkw = 45    # heat conductivity of tube wall (W/m/K)
        self.L = 2       # tube length (m)
        self.Ls = 27     # tube sheet thickness (mm)
        self.rf = 0.025  # tube roughness (mm)
        self.Nt = 86     # total number of tubes in tube bundle
        self.Layout = 1  # tube layout (1:triangular, 2:in-line square, 3:rotated square)
        self.Pt = 1.25 * self.Do  # tube pitch (mm)
        self.Nss = 0     # number of pairs of sealing strips
        self.Npass = 4   # number of passes
        self.Rdt = 0.00036  # fouling resistance (m^2*K/W)
        self.Rds = 0.00018  # fouling resistance (m^2*K/W)
        self.Ds = 305    # shell inside diameter (mm)
        self.Dotl = 294  # shell outside tube limit (mm)
        self.Dsb = 4.45  # shell-baffle clearance (mm)
        self.Lbin = 165  # inlet baffle spacing (mm)
        self.Lbout = 165 # outlet baffle spacing (mm)
        self.Lbc = 450   # central baffle spacing (mm)
        self.Lc = 0.25 * self.Ds  # baffle cut (mm)
        
        # Specification of key variables
        self.Ti1 = 298   # inlet temperatures for tube-side (K)
        self.Ti2 = 303   # outlet temperatures for tube-side (K)
        self.Ts1 = 353   # Shell-side (hot stream) inlet temperature (K)
        self.Ws = 4.8    # Shell-side (hot stream) mass flow rate (kg/sec)
        self.fsT = 1     # state of tube-side fluid (1:liquid, 2:vapor)
        self.fsS = 1     # state of shell-side fluid (1:liquid, 2:vapor)
        self.ptype = 8   # problem type
        self.Ts2 = 338   # guess outlet temperature of shell-side fluid
        self.Wi = 5      # guess tube-side flow rate

# 사용 예시:
# hx_data = HeatExchangerData()
# print(f"Tube Outer Diameter: {hx_data.Do} mm")
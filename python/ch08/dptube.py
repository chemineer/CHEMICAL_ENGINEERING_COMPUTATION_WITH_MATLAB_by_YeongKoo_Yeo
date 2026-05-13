import math

def dptube(D, L, rf, v, rho, Visc, Phi, Npass):
    """
    Calculate tube-side pressure drop (Python version of dptube.m)
    input:
     D: tube inside diameter (mm)
     L: tube length (m)
     rf: tube roughness (mm)
     v: tube-side flow rate (m/s)
     rho: fluid density (kg/m^3)
     Visc: fluid viscosity (Ns/m^2)
     Phi: viscosity correction factor
     Npass: number of tube passes in bundle
    """
    g = 9.81
    rfD = rf / D
    Dm = 1e-3 * D  # mm를 m로 변환
    Nre = Dm * v * rho / Visc
    
    # Friction factor calculation
    if Nre <= 2100:  # Laminar flow
        f = 16 / Nre
    elif Nre <= 4000:  # Zigrang & Sylvester friction factor correlation
        t1 = rfD / 3.7
        t2 = 5.02 / Nre
        # ftm calculation
        ftm = math.log10(t1 - t2 * math.log10(t1 + 13 / Nre))
        f = 1 / (4 * math.log10(t1) + t2 * ftm)**2
    else:  # Round friction factor correlation
        f = 1 / (3.6 * math.log10(Nre / (0.135 * (Nre * rfD + 6.5))))**2
        
    # Pressure drop calculations
    pD = 2 * f * rho * (v**2) * L / Dm
    dp1 = pD / Phi
    dp2 = 4 * rho * L * (v**2) / (2 * g)
    
    DPtube = (dp1 + dp2) * Npass
    
    return DPtube
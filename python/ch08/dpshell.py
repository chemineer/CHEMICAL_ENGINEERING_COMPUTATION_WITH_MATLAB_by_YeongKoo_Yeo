import numpy as np

def dpshell(D, Ds, L, Lbc, Lbin, Lbout, Lc, Dotl, Dsb, Pt, Nt, Nss, w, rho, Visc, Phi, Layout):
    """
    Calculate shell-side pressure drop (Python version of dpshell.m)
    """
    
    # Correlational coefficients for tube arrangement
    # Layouts: 0: triangular, 1: in-line square, 2: rotated square (Python 0-indexing)
    # Original Matlab: 1: triangular, 2: in-line square, 3: rotated square
    layout_idx = Layout - 1
    
    b1 = np.array([
        [0.372, 0.486, 4.570, 45.100, 48.000],
        [0.391, 0.0815, 6.090, 32.100, 35.000],
        [0.303, 0.333, 3.500, 26.200, 32.000]
    ])
    
    b2 = -np.array([
        [0.123, 0.152, 0.476, 0.973, 1.000],
        [0.148, -0.022, 0.602, 0.963, 1.000],
        [0.126, 0.136, 0.476, 0.913, 1.000]
    ])
    
    b3 = np.array([7.00, 6.30, 6.59])
    b4 = np.array([0.500, 0.378, 0.520])

    # Tube layout geometry
    if Layout == 1:
        Pp = 0.866 * Pt
        Pn = Pt / 2
        Pd = Pt
    elif Layout == 2:
        Pp = Pt
        Pn = Pt
        Pd = Pn
    else: # Layout 3 or otherwise
        Pp = 0.7071 * Pt
        Pn = Pp
        Pd = Pn

    # Friction factor for ideal tube-bank
    Nc = Ds * (1 - 2 * Lc / Ds) / Pp
    Ptd = Pt / D
    Sm = Lbc * (Ds - Dotl + (Pt - D) * (Dotl - D) / Pd)
    
    # Shell-side Reynolds number
    Nres = D * w * 1e3 / (Visc * Sm)
    
    # Determine interval J and coefficient c1
    if Nres >= 1e4:
        J = 0 # Matlab index 1
        c1 = 3.7
    elif Nres >= 1e3:
        J = 1 # Matlab index 2
        c1 = 3.7
    elif Nres >= 100:
        J = 2 # Matlab index 3
        c1 = 3.7
    elif Nres >= 10:
        J = 3 # Matlab index 4
        c1 = 4.5
    else:
        J = 4 # Matlab index 5
        c1 = 4.5

    b = b3[layout_idx] / (1 + 0.14 * (Nres**b4[layout_idx]))
    Fj = b1[layout_idx, J] * ((1.33 / Ptd)**b) * (Nres**b2[layout_idx, J])
    
    # Pressure drop for ideal tube bank
    Dpbi = 2 * Fj * Nc * (w / (Sm * 1e-6))**2 / (rho * Phi)

    # Pressure drop for ideal window section
    bdm = (Ds - 2 * Lc) / Dotl
    # Ensure bdm is within [-1, 1] for acos
    bdm = np.clip(bdm, -1.0, 1.0)
    bd1 = np.arccos(bdm)
    bd2 = 1 - 2 * Lc / Ds
    
    Ftc = (np.pi + 2 * bdm * np.sin(bd1) - 2 * bd1) / np.pi
    
    # Ensure value inside arccos is clipped
    bd2_clipped = np.clip(bd2, -1.0, 1.0)
    Sw = (Ds**2 / 4) * (np.arccos(bd2_clipped) - bd2 * np.sqrt(1 - bd2**2)) - (Nt / 8) * (1 - Ftc) * np.pi * D**2
    
    Dw = 4 * Sw / (1.5708 * Nt * (1 - Ftc) * D + 2 * Ds * bd2)
    Ncw = 0.81 * Lc / Pp
    Gw = 1e6 * w / np.sqrt(Sm * Sw)

    if Nres >= 100:
        Dpw = 0.5 * (Gw**2) * (2 + 0.6 * Ncw) / rho
    else:
        Dpw = (2.6e4 * Visc * (Ncw / (Pt - D) + Lbc / (Dw**2)) + Gw) * Gw / rho

    # Correction factor for baffle leakage effects
    Stb = 0.6223 * D * (1 + Ftc) * Nt
    Ssb = Ds * Dsb * 0.5 * (np.pi - np.arccos(bd2_clipped))
    
    R1 = (Stb + Ssb) / Sm
    R2 = Ssb / (Ssb + Stb)
    Pv = -0.15 * (1 + R2) + 0.8
    Rl = np.exp(-1.33 * (1 + R2) * (R1**Pv))

    # Correction factor for bundle bypassing
    Fbp = (Ds - Dotl) * Lbc / Sm
    Nsc = Nss / Nc
    
    if Nsc >= 0.5:
        Rb = 1.0
    else:
        if Nss == 0:
            c2 = 0.0
        else:
            c2 = (2 * Nsc)**0.3333
        Rb = np.exp(-c1 * c2 * Fbp)

    # Correction factor due to unequal baffle spacing
    Nb = 1e3 * L / Lbc + 1
    if Nres >= 100:
        An = 0.2
    else:
        An = 1.0
        
    Rs = (Lbin / Lbc)**(2 - An) + (Lbout / Lbc)**(2 - An)

    # Total shell-side pressure drop
    DPs = Rb * Dpbi * ((Nb - 1) * Rl + (1 + Ncw / Nc) * Rs) + Nb * Rl * Dpw
    
    return DPs
import numpy as np

def htcshell(D, Ds, L, Lbc, Lbin, Lbout, Lc, Dotl, Dsb, Pt, Nt, Nss, w, Visc, Cp, Xk, Phi, Layout):
    """
    Calculate shell-side heat transfer coefficient (Python version of htcshell.m)
    """
    # Correlational coefficients for tube arrangement
    # Layouts: 1:triangular, 2:in-line square, 3:rotated square
    layout_idx = Layout - 1
    
    a1 = np.array([
        [0.321, 0.321, 0.593, 1.360, 1.400],
        [0.370, 0.107, 0.408, 0.900, 0.970],
        [0.370, 0.370, 0.730, 0.498, 1.550]
    ])
    
    a2 = -np.array([
        [0.388, 0.388, 0.477, 0.657, 0.667],
        [0.395, 0.266, 0.460, 0.631, 0.667],
        [0.396, 0.396, 0.500, 0.656, 0.667]
    ])
    
    a3 = np.array([1.450, 1.187, 1.930])
    a4 = np.array([0.519, 0.370, 0.500])

    # Tube layout geometry
    if Layout == 1:
        Pp, Pn, Pd = 0.866 * Pt, Pt / 2, Pt
    elif Layout == 2:
        Pp, Pn, Pd = Pt, Pt, Pt
    else: # Layout 3
        Pp, Pn, Pd = 0.7071 * Pt, Pp, Pp

    Nc = Ds * (1 - 2 * Lc / Ds) / Pp
    Sm = Lbc * (Ds - Dotl + (Pt - D) * (Dotl - D) / Pd)
    Nres = D * w * 1e3 / (Visc * Sm)
    Pr = Cp * Visc / Xk
    Ptd = Pt / D

    # Heat transfer coefficients for ideal tube bank
    if Nres >= 1e4: J, c1 = 0, 1.25
    elif Nres >= 1e3: J, c1 = 1, 1.25
    elif Nres >= 100: J, c1 = 2, 1.25
    elif Nres >= 10: J, c1 = 3, 1.35
    else: J, c1 = 4, 1.35

    a = a3[layout_idx] / (1 + 0.14 * (Nres**a4[layout_idx]))
    Hj = a1[layout_idx, J] * ((1.33 / Ptd)**a) * (Nres**a2[layout_idx, J])
    Hsi = Hj * Cp * Phi * w * (Pr**(-0.6667)) / (Sm * 1e-6)

    # Correction factor for baffle configuration (Phic)
    adm = (Ds - 2 * Lc) / Dotl
    adm = np.clip(adm, -1.0, 1.0)
    adm1 = np.arccos(adm)
    Ftc = (np.pi + 2 * adm * np.sin(adm1) - 2 * adm1) / np.pi
    Phic = Ftc + 0.54 * ((1 - Ftc)**0.345)

    # Correction factor for baffle leakage (Phil)
    Stb = 0.6223 * D * (1 + Ftc) * Nt
    Ssb = Ds * Dsb * 0.5 * (np.pi - np.arccos(np.clip(1 - 2 * Lc / Ds, -1.0, 1.0)))
    R1 = (Stb + Ssb) / Sm
    R2 = Ssb / (Ssb + Stb)
    Phil = 0.44 * (1 - R2) + (1 - 0.44 * (1 - R2)) * np.exp(-2.2 * R1)

    # Correction factor for bundle bypassing (Phib)
    Fbp = (Ds - Dotl) * Lbc / Sm
    Nsc = Nss / Nc
    if Nsc >= 0.5:
        Phib = 1.0
    else:
        c2 = (2 * Nsc)**0.3333 if Nss != 0 else 0
        Phib = np.exp(-c1 * c2 * Fbp)

    Nb = 1e3 * L / Lbc + 1

    # Correction factor for adverse temperature gradient (Phir)
    Phir = 1.0
    if Nres < 100:
        Ncw = 0.8 * Lc / Pp
        Phs = 1.51 / ((Nc + Ncw) * (Nb + 1))**0.18
        if Nres <= 20:
            Phir = Phs
        elif Nres <= 100:
            Phir = Phs - (1 - Phs) * (0.25 - 0.0125 * Nres)
        if Phir <= 0.4:
            Phir = Phs

    # Correction factor for unequal baffle spacing (Phis)
    An = 0.6 if Nres >= 100 else 0.333
    Phis = (Nb - 1 + (Lbin / Lbc)**(1 - An) + (Lbout / Lbc)**(1 - An)) / (Nb - 1 + (Lbin + Lbout) / Lbc)

    Hshell = Hsi * Phic * Phil * Phib * Phir * Phis
    return Hshell
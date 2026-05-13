import numpy as np

def exbco(V, z, Ca0, Fa0, Cp0, Cpc, Ua, m, T1, T2, k1, K2, E, R, dH, hx):
    k = k1 * np.exp(E * (1/T1 - 1/z[2]) / R)
    Kc = K2 * np.exp(dH * (1/T2 - 1/z[2]) / R)
    ra = -k * Ca0 * (1 - (1 + 1/Kc) * z[1])
    dz = np.zeros(3)
    if hx == 'co': dz[0] = Ua * (z[2] - z[0]) / (m * Cpc)
    elif hx == 'cn': dz[0] = -Ua * (z[2] - z[0]) / (m * Cpc)
    else: dz[0] = 0
    dz[1] = -ra / Fa0
    dz[2] = (ra * dH - Ua * (z[2] - z[0])) / (Fa0 * Cp0)
    return dz
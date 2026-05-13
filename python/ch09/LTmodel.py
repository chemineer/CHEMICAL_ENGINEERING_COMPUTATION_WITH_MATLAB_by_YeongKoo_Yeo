import numpy as np

def LTmodel(t, z, ht):
    # z: [h1, h2, T1, T2]
    h1, h2, T1, T2 = z[0], z[1], z[2], z[3]
    
    # Flow rates
    F1b = ht.c1 * np.sqrt(h1 - h2)
    F2 = ht.c2 * np.sqrt(h2)
    dh = max(ht.H, h1) - max(ht.H, h2)
    F1t = ht.c1 * np.sqrt(dh)
    
    # Heat input
    qi1 = ht.Q1 / (ht.rCp * ht.A1 * h1)
    qi2 = ht.Q2 / (ht.rCp * ht.A2 * h2)
    
    # Differential equations
    dz = np.array([
        (ht.F0 - F1t - F1b) / ht.A1,
        (F1t + F1b - F2) / ht.A2,
        (ht.F0 / ht.A1 / h1) * (ht.T0 - T1) + qi1,
        (F1t + F1b) * (T1 - T2) + qi2
    ])
    return dz
import numpy as np

def cstrmult(C, ka, kc, Ca0, Cb0, v0, V):
    # C(0)=Ca, C(1)=Cb, C(2)=Cc, C(3)=Cd (MATLAB의 C(3)=Cb는 원본 오타로 판단하여 Cc로 반영)
    ra = -ka * C[0] * C[1]**2 - 2 * kc * C[0]**2 * C[2]**3 / 3
    rb = -2 * ka * C[0] * C[1]**2
    rc = ka * C[0] * C[1]**2 - kc * C[0]**2 * C[2]**3
    rd = kc * C[0]**2 * C[2]**3 / 3
    
    fC = np.array([
        v0 * Ca0 - v0 * C[0] + ra * V,
        v0 * Cb0 - v0 * C[1] + rb * V,
        -v0 * C[2] + rc * V,
        -v0 * C[3] + rd * V
    ])
    return fC
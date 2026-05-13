import numpy as np
from scipy.optimize import fsolve

def compBCD(Tr, b, c, d, ind):
    # ind는 0 또는 1 (MATLAB의 1, 2에 대응)
    B = b[ind, 0] - b[ind, 1]/Tr - b[ind, 2]/Tr**2 - b[ind, 3]/Tr**3
    C = c[ind, 0] - c[ind, 1]/Tr + c[ind, 2]/Tr**3
    D = d[ind, 0] + d[ind, 1]/Tr
    return B, C, D

def BWReq(Vr, B, C, D, c_ind, Tr, Pr, ind, beta, gam):
    # c_ind는 c[ind, 4]를 의미
    term = c_ind / (Tr**3 * Vr**2) * (beta + gam / Vr**2) * np.exp(-gam / Vr**2)
    return 1 + B/Vr + C/Vr**2 + D/Vr**5 + term - Pr * Vr / Tr

def zLK(T, Tc, P, Pc, w):
    """
    Lee-Kesler 방정식을 사용하여 압축 인자(Z)를 계산합니다.
    """
    b = np.array([[0.1181193, 0.2657280, 0.1547900, 0.0303230],
                  [0.2026579, 0.3315110, 0.0276550, 0.2034880]])
    c = np.array([[0.0236744, 0.0186984, 0.0000000, 0.0427240],
                  [0.0313385, 0.0503618, 0.0169010, 0.0415770]])
    d = 1e-4 * np.array([[0.155488, 0.623689],
                         [0.487360, 0.0740336]])
    beta = np.array([0.653920, 1.22600])
    gam = np.array([0.060167, 0.03754])
    wr = 0.3978
    
    Pr = P / Pc
    Tr = T / Tc
    
    # 1. 단순 유체 (Simple fluid, ind=0)
    ind = 0
    B, C, D = compBCD(Tr, b, c, d, ind)
    Vr0_init = Tr / Pr
    Vr0 = fsolve(BWReq, Vr0_init, args=(B, C, D, c[ind, 3], Tr, Pr, ind, beta[ind], gam[ind]))[0]
    Z0 = Pr * Vr0 / Tr
    
    # 2. 기준 유체 (Reference fluid, ind=1)
    ind = 1
    B, C, D = compBCD(Tr, b, c, d, ind)
    # Vr 초기값은 앞선 계산값 사용
    Vr = fsolve(BWReq, Vr0, args=(B, C, D, c[ind, 3], Tr, Pr, ind, beta[ind], gam[ind]))[0]
    Zr = Pr * Vr / Tr
    
    # 3. 결과 계산
    Z1 = 1/wr * (Zr - Z0)
    Z = Z0 + w * Z1
    
    print(f'Compressibility factor Z = {Z:.6f}')
    return Z
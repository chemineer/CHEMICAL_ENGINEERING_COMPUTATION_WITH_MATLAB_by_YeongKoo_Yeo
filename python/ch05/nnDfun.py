import numpy as np
from scipy.optimize import fsolve

def nnDfun(D, W, rho, dP, rf, L, k, n, K):
    eD = rf / D  # 거칠기
    v = 4 * W / (np.pi * rho * D**2)  # 유속 (m/s)
    
    # 비뉴턴 유체 레이놀즈 수
    Nre = (D**n * v**(2-n) * rho) / (8**(n-1) * k * ((3*n+1)/(4*n))**n)
    
    if Nre < 2100:
        f = 16 / Nre
    else:
        f0 = 16 / Nre
        # Colebrook형 비선형 방정식 해 구하기
        fe = lambda x: 4 * np.log10(Nre * x**(1 - n/2)) / n**0.75 - 0.4 * n**(-1.2) - 1 / np.sqrt(x)
        f = fsolve(fe, f0)[0]
    
    Le = K * D / (4 * f)  # 상당 길이
    # 관경 D를 구하기 위한 잔차 식 반환
    fD = D - ((2 * f * (L + Le)) / (rho * dP) * (4 * W / np.pi)**2)**0.2
    return fD
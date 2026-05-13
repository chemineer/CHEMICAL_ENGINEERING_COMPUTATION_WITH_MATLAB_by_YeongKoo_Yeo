import numpy as np
from scipy.integrate import quad
from deptfun import deptfun

def delHS(state, eos, T1, P1, T2, P2, A, B, C, D, Tc, Pc, w):
    """
    상태 변화에 따른 순수 유체의 엔탈피(H) 및 엔트로피(S) 변화 계산
    """
    # 1. 각 상태에서의 잔류 특성(Residual properties) 계산
    Z1, V1, dH1, dS1 = deptfun(state, eos, T1, P1, Tc, Pc, w)
    Z2, V2, dH2, dS2 = deptfun(state, eos, T2, P2, Tc, Pc, w)
    
    # 2. 이상 기체 상태에서의 변화 계산
    R = 8.314
    # Cp = A + B*T + C*T^2 + D*T^3
    fH = lambda T: A + B*T + C*(T**2) + D*(T**3)
    fS = lambda T: (A/T) + B + C*T + D*(T**2)
    
    dHi, _ = quad(fH, T1, T2)
    dSi_temp, _ = quad(fS, T1, T2)
    dSi = dSi_temp - R * np.log(P2 / P1)
    
    # 3. 최종 H, S 변화량 계산
    dH = dH2 + dHi - dH1
    dS = dS2 + dSi - dS1
    
    return dH, dS

# 예시 호출
# [dH, dS] = delHS('V', 'PR', 300, 1, 400, 2, A, B, C, D, Tc, Pc, w)
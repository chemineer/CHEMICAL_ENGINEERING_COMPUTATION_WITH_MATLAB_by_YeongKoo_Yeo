import numpy as np

def adfun(V, X, pf):
    """
    아세톤 분해 반응의 반응기 설계 방정식 (상미분 방정식)
    
    Parameters:
    V  : 부피 (독립 변수, ode solver의 t 역할)
    X  : [FA, FB, FC, T] 상태 변수 배열
         X[0]: FA, X[1]: FB, X[2]: FC, X[3]: T
    pf : [P, FN2] 시스템 매개변수 배열
    """
    P = pf[0]
    FN2 = pf[1]
    T = X[3]
    
    # 아세톤 농도 계산 (CA)
    # 총 몰수 = FA + FB + FC + FN2
    total_mols = X[0] + X[1] + X[2] + FN2
    CA = 1000 * (X[0] / total_mols) * P / (8.31 * T)
    
    # 속도 상수 및 엔탈피 변화
    k = np.exp(34.34 - 34222 / T)
    dH = 80770 + 6.8 * (T - 298) - 5.75e-3 * (T**2 - 298**2) - 1.27e-6 * (T**3 - 298**3)
    
    # 비열(Cp) 계산
    CpA = 26.2 + 0.183 * T - 45.86e-6 * T**2
    CpB = 20.04 + 0.0945 * T - 30.95e-6 * T**2
    CpC = 13.39 + 0.077 * T - 18.91e-6 * T**2
    CpN2 = 6.25 + 0.00878 * T - 2.1e-8 * T**2
    
    # 반응 속도
    rA = -k * CA
    
    # 상미분 방정식 정의 (dX/dV)
    # dX[0]: dFA/dV, dX[1]: dFB/dV, dX[2]: dFC/dV, dX[3]: dT/dV
    dFA = rA
    dFB = -rA
    dFC = -rA
    dT = -rA * (-dH) / (X[0] * CpA + X[1] * CpB + X[2] * CpC + FN2 * CpN2)
    
    dX = np.array([dFA, dFB, dFC, dT])
    
    return dX
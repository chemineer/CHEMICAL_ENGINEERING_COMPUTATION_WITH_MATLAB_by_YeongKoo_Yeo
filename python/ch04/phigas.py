import numpy as np
from scipy.optimize import fsolve

def phigas(state, eos, T, P, Tc, Pc, w):
    """
    Virial 및 Cubic EOS를 사용하여 압축 인자와 퓨개시티 계수를 추정합니다.
    
    입력:
    state: 상태 ('L': 액체, 'V': 기체)
    eos: 상태 방정식 ('VR', 'VDW', 'RK', 'SRK', 'PR')
    T, P: 온도(K) 및 압력(bar)
    Tc, Pc: 임계 온도(K) 및 임계 압력(bar)
    w: 편심 인자 (acentric factor)
    
    출력:
    phig: 퓨개시티 계수
    f: 퓨개시티 (bar)
    """
    
    # 기본 상수 및 변수 설정
    Tr = T / Tc
    Pr = P / Pc
    R = 83.14  # cm^3*bar/mol/K
    eos = eos.upper()
    state = state.upper()
    
    # EOS별 매개변수 설정
    if eos == 'VDW':
        al, sm, ep, om, ps, kappa = 1, 0, 0, 0.125, 0.42188, 0
    elif eos == 'RK':
        al = 1.0 / np.sqrt(Tr)
        sm, ep, om, ps = 1, 0, 0.08664, 0.42748
    elif eos == 'SRK':
        kappa = 0.480 + 1.574 * w - 0.176 * w**2
        al = (1 + kappa * (1 - np.sqrt(Tr)))**2
        sm, ep, om, ps = 1, 0, 0.08664, 0.42748
    else:  # PR 또는 VR (기본값 PR 세팅)
        kappa = 0.37464 + 1.54226 * w - 0.26992 * w**2
        al = (1 + kappa * (1 - np.sqrt(Tr)))**2
        sm, ep, om = 1 + np.sqrt(2), 1 - np.sqrt(2), 0.0778
        ps = 0.45724

    # 압축 인자 (Z) 계산
    beta = om * Pr / Tr
    q = ps * al / (om * Tr)
    
    if eos == 'VR':  # Virial EOS (기체상만 해당)
        B0 = 0.083 - 0.422 / (Tr**1.6)
        B1 = 0.139 - 0.172 / (Tr**4.2)
        B = R * Tc * (B0 + w * B1) / Pc
        Z = 1 + B * P / (R * T)
    else:
        # Cubic EOS 해를 구하기 위한 함수 정의
        def fV(Z):
            return 1 + beta - q * beta * (Z - beta) / ((Z + ep * beta) * (Z + sm * beta)) - Z
        
        def fL(Z):
            return beta + (Z + ep * beta) * (Z + sm * beta) * (1 + beta - Z) / (q * beta) - Z
        
        if state == 'V':
            Z = fsolve(fV, 1.0)[0]
        elif state == 'L':
            Z = fsolve(fL, beta if beta > 0 else 1e-5)[0]
        else:
            Z = 1.0

    # 몰 부피 계산
    V = Z * R * T / P  # cm^3/mol
    
    # 퓨개시티 계수 (phig) 계산
    a = ps * al * (R**2) * (Tc**2) / Pc
    b = om * R * Tc / Pc
    qi = a / (b * R * T)
    Bd = b * P / (R * T)
    
    if eos == 'VR':
        B0 = 0.083 - 0.422 / (Tr**1.6)
        B1 = 0.139 - 0.172 / (Tr**4.2)
        phig = np.exp(Pr * (B0 + w * B1) / Tr)
    elif eos == 'VDW':
        phig = np.exp(Z - 1 - np.log(Z * (1 - b/V)) - a / (R * T * V))
    elif eos == 'RK':
        phig = np.exp(Z - 1 - np.log(Z * (1 - b/V)) - (a / (b * R * T)) * np.log(1 + b/V))
    else:  # SRK 또는 PR
        phig = np.exp(Z - 1 - np.log(Z - Bd) - (qi / (sm - ep)) * np.log((Z + sm * Bd) / (Z + ep * Bd)))
        
    f = phig * P  # 퓨개시티 (bar)
    
    return phig, f

# 사용 예시 (예: Water, PR, Vapor)
# phig_val, f_val = phigas('V', 'PR', 500, 50, 647.1, 220.55, 0.344)
# print(f"Fugacity Coefficient: {phig_val:.4f}, Fugacity: {f_val:.4f} bar")
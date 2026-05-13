import numpy as np
from scipy.optimize import fsolve

def deptfun(state, eos, T, P, Tc, Pc, w):
    """
    비리얼 및 입방 상태 방정식을 이용한 잔류 성질(Departure Functions) 계산
    
    Returns:
    Z: 압축 인자
    V: 몰 부피 (cm^3/mol)
    dH: 엔탈피 잔류 (J/mol)
    dS: 엔트로피 잔류 (J/mol/K)
    """
    # 기본 상수 및 환산 변수 설정
    Tr = T / Tc
    Pr = P / Pc
    R_bar = 83.14  # cm^3*bar/mol/K
    eos = eos.upper()
    state = state.upper()
    
    # EOS별 파라미터 설정
    if eos == 'VDW':
        al, sm, ep, om, ps, kappa = 1, 0, 0, 0.125, 0.42188, 0
    elif eos == 'RK':
        al = 1.0 / np.sqrt(Tr)
        sm, ep, om, ps = 1, 0, 0.08664, 0.42748
        kappa = 0 # RK는 kappa를 사용하지 않음
    elif eos == 'SRK':
        kappa = 0.480 + 1.574*w - 0.176*w**2
        al = (1 + kappa * (1 - np.sqrt(Tr)))**2
        sm, ep, om, ps = 1, 0, 0.08664, 0.42748
    else:  # PR 또는 VR
        kappa = 0.37464 + 1.54226*w - 0.26992*w**2
        al = (1 + kappa * (1 - np.sqrt(Tr)))**2
        sm, ep, om, ps = 1 + np.sqrt(2), 1 - np.sqrt(2), 0.0778, 0.45724

    # 압축 인자(Z) 계산
    beta = om * Pr / Tr
    q = ps * al / (om * Tr)
    
    if eos == 'VR': # Virial EOS (Vapor phase)
        B0 = 0.083 - 0.422 / (Tr**1.6)
        B1 = 0.139 - 0.172 / (Tr**4.2)
        B = R_bar * Tc * (B0 + w * B1) / Pc
        Z = 1 + B * P / (R_bar * T)
    else: # Cubic EOS (VDW, RK, SRK, PR)
        # fzero를 fsolve로 대체
        if state == 'V':
            fV = lambda Z: 1 + beta - q * beta * (Z - beta) / ((Z + ep * beta) * (Z + sm * beta)) - Z
            Z = fsolve(fV, 1.0)[0]
        else: # Liquid
            fL = lambda Z: beta + (Z + ep * beta) * (Z + sm * beta) * (1 + beta - Z) / (q * beta) - Z
            Z = fsolve(fL, beta)[0]

    V = Z * R_bar * T / P # cm^3/mol
    
    # Departure function 계산을 위한 보조 변수
    a = ps * al * (R_bar**2) * (Tc**2) / Pc
    b = om * R_bar * Tc / Pc
    Ad = a * P / ((R_bar**2) * (T**2))
    Bd = b * P / (R_bar * T)
    
    # EOS별 dH, dS 로직
    if eos == 'VR':
        dH = -Pr * (1.0972 / Tr**2.6 - 0.083 / Tr + w * (0.8944 / Tr**5.2 - 0.139 / Tr)) * R_bar * T
        dS = -Pr * (0.675 / Tr**2.6 + w * 0.722 / Tr**5.2) * R_bar
    elif eos == 'RK':
        dH = (Z - 1 - 1.5 * ps / (om * Tr**1.5) * np.log(1 + b / V)) * R_bar * T
        dS = (np.log(Z - Bd) - 0.5 * ps / (om * Tr**1.5) * np.log(1 + b / V)) * R_bar
    elif eos == 'VDW':
        dH = (Z - 1 - 3.375 * Bd / Tr / Z) * R_bar * T
        dS = R_bar * np.log(Z - Bd)
    elif eos == 'SRK':
        term1 = -kappa * np.sqrt(Tr) / (1 + kappa * (1 - np.sqrt(Tr)))
        dH = (Z - 1 + (term1 - 1) * (Ad / Bd) * np.log(1 + Bd / Z)) * R_bar * T
        dS = (np.log(Z - Bd) + term1 * (Ad / Bd) * np.log(1 + Bd / Z)) * R_bar
    else: # PR
        sqrt8 = np.sqrt(8)
        term_pr = (Ad / (Bd * sqrt8))
        term_k = (kappa * np.sqrt(Tr) / np.sqrt(al))
        log_term = np.log((Z + sm * Bd) / (Z + ep * Bd))
        dH = (Z - 1 - term_pr * (1 + term_k) * log_term) * R_bar * T
        dS = (np.log(Z - Bd) - term_pr * term_k * log_term) * R_bar

    # 단위 변환: R_bar(83.14)에서 R(8.314) 단위계로 조정 (1/10배)
    dH = dH / 10.0 # J/mol
    dS = dS / 10.0 # J/mol/K
    
    return Z, V, dH, dS
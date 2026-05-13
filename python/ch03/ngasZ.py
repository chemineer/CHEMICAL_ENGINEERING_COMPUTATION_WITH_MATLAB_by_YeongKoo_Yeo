import numpy as np

def ngasZ(T, P, Sg):
    """
    Estimates the compressibility factor Z of natural gases
    T: temperature (F) (scalar or vector)
    P: pressure (psia)
    Sg: specific gravity of the natural gas
    """
    
    # 단위 변환: psia를 1000psia 단위로, 화씨를 랭킨(Rankine)으로 변환
    P_val = P / 1000.0
    T_val = T + 460.0
    
    # 상수 정의
    A1, A2, A3, A4, A5, A6 = 0.001946, -0.027635, 0.136315, -0.23849, 0.105168, 3.44e8
    
    # 계산식
    F1 = P_val * (0.251 * Sg - 0.15) - 0.202 * Sg + 1.106
    
    den = 1 + (A6 * P_val * (10.0 ** (1.785 * Sg))) / (T_val ** 3.825)
    
    F2 = 1.4 * np.exp(-0.0054 * (T_val - 460.0))
    
    F3 = A1 * (P_val ** 5) + A2 * (P_val ** 4) + A3 * (P_val ** 3) + A4 * (P_val ** 2) + A5 * P_val
    
    F4 = (0.154 - 0.152 * Sg) * (P_val ** (3.18 * Sg - 1)) * np.exp(-0.5 * P_val) - 0.02
    
    F5 = 0.35 * (0.6 - Sg) * np.exp(-1.039 * (P_val - 1.8) ** 2)
    
    # 결과 계산
    nz = F1 * (1.0 / den + F2 * F3) + F4 + F5
    
    return nz
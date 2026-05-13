import numpy as np
from scipy.integrate import quad

def delH(C, T1, T2):
    """
    순수 물질의 엔탈피 변화(Delta H)와 평균 열용량을 계산합니다.
    
    Parameters:
    C  : Cp 관계식 계수 [A, B, C, D]
    T1 : 하한 온도 (K)
    T2 : 상한 온도 (K)
    
    Returns:
    Q  : 엔탈피 변화 (J)
    mc : 평균 열용량 비율 ((Cp)h / R)
    """
    R = 8.314  # J/(mol-K)
    
    # Cp 관계식 정의 (C=[A, B, C, D])
    # fH(T) = A + B*T + C*T^2 + D*T^-2
    def fH(T):
        return C[0] + C[1]*T + C[2]*(T**2) + C[3]*(T**(-2))
    
    # 수치 적분 (MATLAB의 quadl 대체)
    # quad는 (적분값, 오차추정값)을 반환하므로 [0]번 인덱스만 사용
    intCp, error = quad(fH, T1, T2)
    
    Q = R * intCp  # J
    mc = intCp / (T2 - T1)  # (Cp)h/R
    
    return Q, mc
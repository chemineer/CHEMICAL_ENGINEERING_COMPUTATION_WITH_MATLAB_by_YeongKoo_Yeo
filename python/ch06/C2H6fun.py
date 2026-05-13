import numpy as np

def C2H6fun(V, F, P, T, Fs):
    """
    에탄 분해 반응의 몰 수지 방정식
    
    Parameters:
    V  : 반응기 부피 (독립 변수)
    F  : [F1, F2, F3, F4, F5, F6, F7, F8] (몰 유량 배열)
    P  : 압력
    T  : 온도
    Fs : 추가 성분(예: 희석제)의 몰 유량
    """
    # 기체 상수
    Re = 1.987
    Ri = 0.08314
    
    # 반응 속도 상수 (Arrhenius 식)
    k1 = 4.65e13 * np.exp(-65210 / (Re * T))
    k2 = 3.85e11 * np.exp(-65210 / (Re * T))
    k3 = 9.81e8 * np.exp(-36920 / (Re * T))
    k4 = 1.03e12 * np.exp(-41260 / (Re * T))
    k5 = 7.08e13 * np.exp(-60430 / (Re * T))
    
    # 농도 계산
    Ft = np.sum(F) + Fs
    C = F * P / (Ft * Ri * T)
    
    # 반응 속도 정의
    r1 = k1 * C[0]
    r2 = k2 * C[0]**2
    r3 = k3 * C[3]
    r4 = k4 * C[1] * C[2]
    r5 = k5 * C[0] * C[1]
    
    # 몰 수지 식 (dF/dV)
    dF = np.array([
        -r1 - 2*r2 - r5,      # C2H6 (1)
        r1 - r4 - r5,         # C2H4 (2)
        r3 - r4,              # C2H2 (3)
        -r3 + r5,             # C3H6 (4)
        r2,                   # C3H8 (5)
        r4,                   # C4H6 (6)
        r2 + r3 + r5,         # CH4 (7)
        r1                    # H2 (8)
    ])
    
    return dF
import numpy as np

def vpRM(T, Tb, nu, dB, GI, Dp):
    """
    Rarey/Moller 방정식을 이용한 증기압(vapor pressure) 추정
    
    Parameters:
    T  : 온도 (K) (스칼라 또는 넘파이 배열)
    Tb : 일반 끓는점 (K)
    nu : 그룹의 빈도 (벡터)
    dB : delta B 값 (벡터)
    GI : GIij 값 (행렬 또는 배열)
    Dp : D prime 값
    
    Returns:
    Pv : 증기압 (bar)
    """
    # 입력값을 넘파이 배열로 변환
    nu = np.array(nu)
    dB = np.array(dB)
    GI = np.array(GI)
    
    # 합계 계산
    sumdBi = np.sum(nu * dB)
    sumGI = np.sum(GI)
    
    # Bp 계산
    Bp = 9.42208 + sumdBi + sumGI
    
    # 증기압(Pv) 계산 (bar)
    # MATLAB의 (Tb.^1.485)는 배열 연산을 고려하여 작성됨
    Pv = np.exp(Bp * (T - Tb) / (T + 2.65 - (Tb**1.485) / 135) + Dp * np.log(T / Tb))
    
    print(f'Vapor pressure = {Pv} bar')
    return Pv
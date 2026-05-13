import numpy as np

def vpwagner(T, Tc, Pc, C):
    """
    Wagner 방정식을 이용한 증기압(Pv) 추정
    
    Parameters:
    T  : 온도 (K) (스칼라 또는 넘파이 배열)
    Tc : 임계 온도 (K)
    Pc : 임계 압력 (MPa)
    C  : Wagner 방정식 파라미터 벡터 [a, b, c, d]
    
    Returns:
    Pv : 추정된 증기압 (MPa)
    """
    # 환산 온도(Reduced Temperature) 계산
    Tr = T / Tc
    
    # 파라미터 분할
    a, b, c, d = C[0], C[1], C[2], C[3]
    
    # Wagner 방정식 계산
    # Pv = Pc * exp( (a(1-Tr) + b(1-Tr)^1.5 + c(1-Tr)^2.5 + d(1-Tr)^5) / Tr )
    exponent = (a * (1 - Tr) + 
                b * (1 - Tr)**1.5 + 
                c * (1 - Tr)**2.5 + 
                d * (1 - Tr)**5) / Tr
    
    Pv = Pc * np.exp(exponent)
    
    print(f'Vapor pressure = {Pv} MPa')
    return Pv
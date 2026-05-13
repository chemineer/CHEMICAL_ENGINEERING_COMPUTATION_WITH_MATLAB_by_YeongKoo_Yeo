import numpy as np

def vtfun(x, rp, ro, mu, dp):
    """
    입자의 종단 속도(x)를 구하기 위한 비선형 방정식 계산
    
    Parameters:
    x  : 종단 속도 (terminal velocity)
    rp : 입자 밀도 (particle density)
    ro : 유체 밀도 (fluid density)
    mu : 유체 점도 (dynamic viscosity)
    dp : 입자 직경 (particle diameter)
    """
    g = 9.8
    # 레이놀즈 수 계산
    Nre = dp * ro * x / mu
    
    # 항력 계수(Cd) 계산 (레이놀즈 수 범위에 따른 조건문)
    if Nre < 0.1:
        Cd = 24.0 / Nre
    elif Nre < 1000:
        Cd = (1 + 0.14 * Nre**0.7) * 24.0 / Nre
    elif Nre < 3.5e5:
        Cd = 0.44
    else:
        Cd = 0.19 - 8e4 / Nre
        
    # 종단 속도 조건에 따른 비선형 방정식 반환
    fvt = 3 * Cd * ro * x**2 - 4 * g * (rp - ro) * dp
    
    return fvt
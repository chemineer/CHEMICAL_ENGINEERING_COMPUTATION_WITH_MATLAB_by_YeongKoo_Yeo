import numpy as np

def vhwagner(T, Tc, Pc, vL, vV, C):
    """
    Estimation of enthalpy of vaporization using the Wagner equation
    
    Parameters:
    T, Tc : Temperature and critical temperature (K)
    Pc    : Critical pressure (MPa)
    vL    : Molar volume of liquid
    vV    : Molar volume of vapor
    C     : Parameter vector of the Wagner equation [a, b, c, d]
    
    Returns:
    dHv   : Estimated enthalpy of vaporization (J/mol)
    """
    
    # 환산 온도 계산
    Tr = T / Tc 
    
    # 파라미터 추출 (파이썬은 인덱스가 0부터 시작함)
    a = C[0]
    b = C[1]
    c = C[2]
    d = C[3]
    
    # Wagner 식을 이용한 증기압(Pv) 계산
    # Pv = Pc * exp((a*(1-Tr) + b*(1-Tr)^1.5 + c*(1-Tr)^2.5 + d*(1-Tr)^5) / Tr)
    exponent = (a * (1 - Tr) + 
                b * (1 - Tr)**1.5 + 
                c * (1 - Tr)**2.5 + 
                d * (1 - Tr)**5) / Tr
    Pv = Pc * np.exp(exponent)
    
    # wd 계산 (중간 계산식)
    wd = (np.log(Pv / Pc) + a + 
          1.5 * b * (1 - Tr)**0.5 + 
          2.5 * c * (1 - Tr)**1.5 + 
          5 * d * (1 - Tr)**4)
    
    # 증발 엔탈피 계산 (J/mol)
    # dHv = -Pv * (vV - vL) * wd * 1e6
    dHv = -Pv * (vV - vL) * wd * 1e6
    
    return dHv
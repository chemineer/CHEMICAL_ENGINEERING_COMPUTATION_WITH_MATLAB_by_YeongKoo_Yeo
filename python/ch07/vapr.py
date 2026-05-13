import numpy as np
from scipy.optimize import bisect

# 이전 단계에서 정의한 vapdat 데이터를 가져오는 함수 (가정)
def get_vapdat():
    return {
        'Q': 6000, 'lambda': 475, 'Cp': 1, 'Ti': 20, 'Fi': 1000, 'rhol': 1,
        'V': 5500, 'Kv': 40, 'P0': 10, 'Mw': 18, 'R': 0.08206,
        'A': 18.3036, 'B': 3816.44, 'C': -46.13
    }

def vapr(t, z):
    """
    증기 시스템 미분 방정식 함수
    z[0] = mV (증기 질량), z[1] = mL (액체 질량)
    """
    # 1. 데이터 로드
    d = get_vapdat()
    A, B, C = d['A'], d['B'], d['C']
    Q, lam, Cp, Ti, Fi, rhol = d['Q'], d['lambda'], d['Cp'], d['Ti'], d['Fi'], d['rhol']
    V, Kv, P0, Mw, R = d['V'], d['Kv'], d['P0'], d['Mw'], d['R']
    
    mV = z[0]
    mL = z[1]
    Vv = V - mL / rhol
    
    # 2. 이분법을 이용한 온도 T 결정
    # f(T) = Psat(T) - P_ideal_gas(T) = 0이 되는 T를 찾음
    def f_temp(T):
        return np.exp(A - B / (T + C)) - (7.6e5 * mV * R * T) / (Vv * Mw)
    
    Ta = 273.15
    Tb = Ta + 400
    
    # MATLAB의 bisection 루프를 scipy bisect로 대체하거나 직접 구현
    try:
        # f(Ta)*f(Tb) > 0인 경우 해가 없을 수 있음
        if f_temp(Ta) * f_temp(Tb) > 0:
            print("No solution T")
            return [0, 0]
        
        T = bisect(f_temp, Ta, Tb, xtol=1e-2)
    except ValueError:
        print("No solution T in the given range")
        return [0, 0]
    
    # 3. 미분 방정식 정의
    # P 계산 (mmHg -> bar 변환을 위한 750.044 나누기)
    P = np.exp(A - B / (T + C)) / 750.044 # P: bar
    
    # vB: 생성되는 증기 유량, vO: 유출되는 증기 유량
    vB = (Fi * rhol * Cp * Ti + Q) / (lam + Cp * (T - 273.15))
    vO = Kv * np.sqrt(max(0, P * (P - P0))) # sqrt 안의 음수 방지
    
    # dz[0] = dmV/dt, dz[1] = dmL/dt
    dz = np.array([vB - vO, Fi * rhol - vB])
    
    return dz
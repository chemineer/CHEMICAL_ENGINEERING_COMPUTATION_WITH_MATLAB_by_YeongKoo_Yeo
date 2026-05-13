import numpy as np
from scipy.optimize import fsolve
from satsteam import satsteam

def sinevap(evdat):
    """
    단일 효과 증발기(Single-effect evaporator) 계산 함수
    """
    # 상수 및 데이터 준비
    Tc = 647.096  # K
    Pc = 22064000  # Pa
    T0 = 273.15
    xf = evdat['xf']
    xp = evdat['xp']
    mf = evdat['mf']
    Tf = (evdat['Tf'] - 32) / 1.8 + T0
    Ps = evdat['Ps'] * 6894.757
    Pv = evdat['Pv'] * 6894.757
    U = evdat['U']
    cf = 1 / 2326  # J/kg -> Btu/lb
    
    # Saturation temperature 계산을 위한 다항식 계수
    a = np.array([-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502])
    
    def get_x(P):
        # f(x) = 0을 만족하는 x를 찾기 위한 함수
        func = lambda x: (a[0]*x + a[1]*x**1.5 + a[2]*x**3 + a[3]*x**3.5 + 
                          a[4]*x**4 + a[5]*x**7.5 - (1-x)*np.log(P/Pc))
        return fsolve(func, 0.5)[0]
    
    xs = get_x(Ps)
    xv = get_x(Pv)
    Ts_k = Tc * (1 - xs)
    Tv_k = Tc * (1 - xv)
    
    # 엔탈피 계산 (외부 satsteam 함수 호출)
    sts = satsteam(Ts_k - T0)
    stv = satsteam(Tv_k - T0)
    stf = satsteam(Tf - T0)
    
    Hs = sts['hV'] * cf
    hs = sts['hL'] * cf
    Hv = stv['hV'] * cf
    hv = stv['hL'] * cf
    hf = stf['hL'] * cf
    
    # 물질 수지 및 열 수지 계산
    mp = xf * mf / xp
    mev = mf - mp
    q = mev * Hv - mf * hf + mp * hv
    ms = q / (Hs - hs)
    
    # 결과 계산 (단위 변환)
    Ts_f = (Ts_k - T0) * 1.8 + 32
    Tv_f = (Tv_k - T0) * 1.8 + 32
    dT = Ts_f - Tv_f
    A = q / U / dT
    
    # 결과 반환
    res = {
        'ms': ms,
        'A': A
    }
    return res
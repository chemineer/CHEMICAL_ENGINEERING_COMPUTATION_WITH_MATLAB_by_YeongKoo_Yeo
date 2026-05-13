import numpy as np
from scipy.integrate import quad
from crfdata import crfdata

def mArea(i, D1, E1, F1, R1, S1, T1, alpa, r, xf, uf):
    """막 면적 계산을 위한 피적분 함수 (MATLAB의 mArea 대응)"""
    u = -D1 * i + np.sqrt((D1**2) * i**2 + 2 * E1 * i + F1**2)
    # 계산 로직 (MATLAB 원본의 mArea 수식을 기반으로 작성)
    term1 = (1 + (1 - alpa) * u / (alpa * (1 - r))) / (u - r)
    term2 = ((1 - xf) / (1 - xf/(1+i))) * ((uf - E1/D1) / (u - E1/D1))**R1
    term3 = ((uf - alpa + F1) / (u - alpa + F1))**S1 * ((uf - F1) / (u - F1))**T1
    return term1 * term2 * term3

def crf1ex():
    """
    십자 흐름 모델을 이용한 이성분 공급물 막 분리 공정 계산
    MATLAB의 crf1ex 함수를 파이썬으로 구현
    """
    # 1. 데이터 로드 (이전에 정의한 crfdata 함수 호출)
    data = crfdata()
    t = data['t']; Pm = data['Pm']; alpa = data['alpa']
    ph = data['ph']; pl = data['pl']; r = data['r']
    qf = data['qf']; xf = data['xf']; xr = data['xr']
    
    Pa = Pm[0]
    Pb = Pm[1]
    
    # 2. 중간 변수 계산
    D1 = ((1 - alpa) * r + alpa) / 2
    F1 = -((1 - alpa) * r - 1) / 2
    E1 = alpa / 2 - D1 * F1
    R1 = 1.0 / (2 * D1 - 1)
    S1 = (alpa * (D1 - 1) + F1) / ((2 * D1 - 1) * (alpa / 2 - F1))
    T1 = 1.0 / (1 - D1 - E1 / F1)
    
    # 성분비 계산
    i0 = xr / (1 - xr)
    i2 = xf / (1 - xf)
    
    # u 값 계산
    ur = -D1 * i0 + np.sqrt((D1**2) * i0**2 + 2 * E1 * i0 + F1**2)
    uf = -D1 * i2 + np.sqrt((D1**2) * i2**2 + 2 * E1 * i2 + F1**2)
    
    # 3. Stage-cut (theta) 및 투과물 농도 (yp) 계산
    theta = 1 - ((1 - xf) / (1 - xr)) * \
            ((uf - E1 / D1) / (ur - E1 / D1))**R1 * \
            ((uf - alpa + F1) / (ur - alpa + F1))**S1 * \
            ((uf - F1) / (ur - F1))**T1
    
    yp = (xf - (1 - theta) * xr) / theta
    
    # 4. 수치 적분을 이용한 막 면적 계산
    # quad 함수를 사용하여 i0부터 i2까지 mArea를 적분
    Ami, error = quad(mArea, i0, i2, args=(D1, E1, F1, R1, S1, T1, alpa, r, xf, uf))
    
    Am = Ami * qf * t / (ph * Pb)
    rc = theta * yp / xf # 회수율 (recovery ratio)
    
    # 결과 반환 (딕셔너리 형태)
    res = {
        'yp': yp,
        'theta': theta,
        'Am': Am,
        'rc': rc,
        'xr': xr
    }
    
    return res

# 참고: crfdata() 함수가 사전에 정의되어 있어야 합니다.
import numpy as np
from cpmdata import cpmdata

def cpm1ex():
    """
    MATLAB의 cpm1ex 함수를 파이썬으로 구현.
    참고: cpmdata()가 먼저 호출되어 Pm, alpa, r, xr, xf, qf, t, ph, pl 등의 
    변수가 정의되어 있어야 합니다.
    """
    # 외부 데이터 로드 (cpmdata 함수가 별도로 정의되어 있다고 가정)
    # 만약 cpmdata를 직접 불러올 수 없다면, 해당 변수들을 인자로 넘겨받도록 수정해야 합니다.
    data = cpmdata() 
    
    # 데이터 추출 (딕셔너리 혹은 별도 객체에서 가져오기)
    Pm = data['Pm']; alpa = data['alpa']; r = data['r']; xr = data['xr']
    xf = data['xf']; qf = data['qf']; t = data['t']; ph = data['ph']; pl = data['pl']

    Pa = Pm[0]  # MATLAB의 Pm(1) 대응
    a = 1 - alpa
    b = -1 + alpa + 1/r + (alpa - 1) * xr / r
    c = -alpa * xr / r
    
    # 2차 방정식의 근의 공식: permeate mole fraction
    yp = (-b + np.sqrt(b**2 - 4 * a * c)) / (2 * a)
    
    # stage-cut
    theta = (xf - xr) / (yp - xr)
    
    # membrane area
    Am = (theta * qf * yp) / ((Pa / t) * (ph * xr - pl * yp))
    
    # recovery ratio
    rc = qf * theta * yp / (qf * xf)
    
    # xom 계산
    xom = (xf * (1 + r * (alpa - 1) * (1 - xf))) / (xf * (1 - alpa) + alpa)
    
    # 결과 객체 반환 (딕셔너리)
    res = {
        'yp': yp,
        'theta': theta,
        'Am': Am,
        'rc': rc,
        'xr': xr
    }
    
    return res


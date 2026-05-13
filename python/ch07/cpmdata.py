import numpy as np

def cpmdata():
    """
    MATLAB의 cpmdata.m 스크립트를 파이썬 함수로 구현.
    반환값은 각 데이터가 담긴 딕셔너리입니다.
    """
    data = {}
    
    # 데이터 할당
    data['t'] = 0.00254       # membrane thickness (cm)
    data['Pm'] = np.array([50, 5]) * 1e-10  # permeability (cm^3*cm/(sec*cm^2*cmHg))
    data['alpa'] = data['Pm'][0] / data['Pm'][1]
    data['ph'] = 80           # feed-side pressure (cmHg)
    data['pl'] = 20           # permeate-side pressure (cmHg)
    data['r'] = data['pl'] / data['ph']  # pressure ratio
    data['qf'] = 1e4          # feed flow rate (cm^3/sec(STP))
    data['xf'] = 0.5          # feed composition (mole fraction)
    data['xr'] = 0.25         # desired reject composition (mole fraction)
    
    return data

# 사용 예시:
# data = cpmdata()
# print(f"Alpha value: {data['alpa']}")
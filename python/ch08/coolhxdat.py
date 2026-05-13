import numpy as np

def coolhxdat():
    """
    열교환기(Heat Exchanger) 관련 데이터를 정의하는 함수입니다.
    데이터를 딕셔너리 형태로 반환합니다.
    """
    # Physical properties
    data = {
        'Cph': 0.92, 'Cpc': 1.0, 'rhoh': 59.76, 'rhoc': 59.87, 
        'kh': 0.3, 'kc': 0.36, 'kw': 30.0, 'muh': 0.75, 
        'muc': 0.77, 'muwc': 0.77, 'muwh': 0.75,
        
        # Operating conditions
        'mh': 6.2e4, 'ui': 5.0, 'T1': 150.0, 'T2': 120.0, 
        't1': 75.0, 't2': 100.0,
        
        # Shell and tube
        'L': 15.0, 'Do': 0.0625, 'Di': 0.04017, 'Ds': 1.4375, 
        'Pt': 1/12, 'cl': 0.02083, 'Rfi': 0.004, 'Rfo': 0.0
    }
    
    # Parameters and basic properties (계산된 항목들)
    data['gc'] = 32.2
    
    # 로그 평균 온도차 (LMTD)
    dTlm = ((data['T1'] - data['t2']) - (data['T2'] - data['t1'])) / np.log((data['T1'] - data['t2']) / (data['T2'] - data['t1']))
    data['dTlm'] = dTlm
    
    # 보정 계수 (Correction factor) 계산
    R = (data['T1'] - data['T2']) / (data['t2'] - data['t1'])
    S = (data['t2'] - data['t1']) / (data['T1'] - data['t1'])
    
    F12den = (R - 1) * np.log((2 - S*(R + 1 - np.sqrt(R**2 + 1))) / (2 - S*(R + 1 + np.sqrt(R**2 + 1))))
    F12 = np.sqrt(R**2 + 1) * np.log((1 - S) / (1 - R*S)) / F12den
    data['F12'] = F12
    
    # 열 부하 및 냉각수 유량
    Q = data['mh'] * data['Cph'] * (data['T1'] - data['T2'])
    mc = Q / (data['Cpc'] * (data['t2'] - data['t1']))
    data['Q'] = Q
    data['mc'] = mc
    
    # 기하학적 파라미터 계산
    Aci = data['mh'] / (data['rhoh'] * data['ui'] * 3600)
    Nt = np.ceil(4 * Aci / (np.pi * data['Di']**2))
    At = np.pi * data['Di'] * data['L']
    
    data['Nt'] = Nt
    data['At'] = At
    data['dx'] = (data['Do'] - data['Di']) / 2
    data['Dlm'] = (data['Do'] - data['Di']) / np.log(data['Do'] / data['Di'])
    data['De'] = (4 / (np.pi * data['Do'])) * (data['Pt']**2 - np.pi * data['Do']**2 / 4)
    data['Ui'] = 150.0 # 가정값
    
    return data
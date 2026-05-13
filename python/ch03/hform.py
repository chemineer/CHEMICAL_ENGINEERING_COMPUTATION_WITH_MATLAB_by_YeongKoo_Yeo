import numpy as np
from compID import compID

def hform(T, cname):
    """
    Estimation of heat of formation of gases (kcal/gmol) at low temperature
    T: temperature(C) (scalar or list/numpy array)
    cname: name or chemical formula of the compound
    """
    
    # correlation constants (A, B, C)
    # 0.0 values indicate missing components in the original matrix
    cv = np.array([
        [0.0, 0.0, 0.0],    # F2
        [0.0, 0.0, 0.0],    # Cl2
        [-69.6000, -5.2900, 0], [-26.5000, 0.8400, -1.0500], [-93.9000, -0.4300, 0], [-21.9000, -0.6100, 0],
        [-9.3400, -6.2000, 2.3600], [-57.4000, -1.7900, 0], [-31.8000, -3.0400, 1.1900],
        [0.0, 0.0, 0.0],    # H2
        [0.0, 0.0, 0.0],    # N2
        [0.0, 0.0, 0.0],    # O2
        [0.0, 0.0, 0.0],    # C2H4
        [-15.4000, -9.5900, 3.5000], [-16.4000, -14.8000, 6.1300], [-20.0000, -19.1000, 8.1500],
        [23.7000, -15.3000, 6.2700], [16.8000, -19.0000, 7.8400], [25.2000, -17.6000, 8.9800],
        [-19.3000, -14.6000, 7.1800], [16.9000, -16.7000, 7.9000], [-21.6000, -32.0000, 15.8000],
        [29.0000, -10.1000, 4.1300], [-44.9600, -11.9000, 4.9800], [-24.7000, 0.0334, 0],
        [-24.0000, 2.4200, 0]
    ])
    
    # 단위 변환 (1e-3, 1e-6 적용)
    wv = np.column_stack([cv[:, 0], cv[:, 1] * 1e-3, cv[:, 2] * 1e-6])
    
    # compID 함수는 외부에서 정의되어 있어야 합니다.
    ind = compID(cname) 
    
    # 온도 변환 (C -> K)
    T = np.array(T)
    T_k = T + 273.15
    
    # Heat of formation 계산
    if ind != 3:
        hf = wv[ind-1, 0] + wv[ind-1, 1] * T_k + wv[ind-1, 2] * (T_k**2)
    else:
        # ind = 3: sulfur dioxide 특수 케이스
        hf = []
        for t in T_k:
            if 298 <= t <= 717:
                hfv = wv[0, 0] + wv[0, 1] * t + wv[0, 2] * (t**2)
            else:
                # 717K~1500K 범위
                wv_so2 = np.array([-86.9, 0.32 * 1e-3, 0])
                hfv = wv_so2[0] + wv_so2[1] * t + wv_so2[2] * (t**2)
            hf.append(hfv)
        hf = np.array(hf)
        
    return hf
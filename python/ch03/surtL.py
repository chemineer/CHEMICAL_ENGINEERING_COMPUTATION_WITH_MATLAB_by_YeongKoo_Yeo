import numpy as np
from compID import compID

def surtL(T, cname):
    """
    Calculates liquid surface tension ratio (sg = sigma/sigma1)
    T: temperature (C) (scalar or row vector)
    cname: common name or chemical formula of the compound
    """
    
    # 데이터 행렬 (T1(C), Tc(C), r)
    wv = np.array([
        [-200.0, -129.0, 0.8811], [20.0, 144.0, 1.0508], [30.0, 157.6, 1.1768],
        [-193.0, -140.1, 1.1441], [20.0, 31.1, 1.3015], [-93.0, 51.5, 1.0972],
        [-45.0, 132.4, 1.1548], [25.0, 374.2, 0.8105], [18.2, 455.0, 0.9141],
        [-256.0, -240.2, 1.1012], [-203.0, -146.8, 1.2123], [-202.0, -118.5, 1.1933],
        [-120.0, 9.9, 1.2760], [-168.16, -82.6, 1.3941], [-120.0, 32.3, 1.2060],
        [-90.0, 96.7, 1.1982], [20.0, 288.94, 1.2243], [20.0, 318.8, 1.2364],
        [20.0, 426.0, 1.1022], [60.0, 420.0, 1.0725], [10.0, 124.9, 1.3201],
        [20.0, 280.3, 1.4246], [40.0, 152.0, 1.2055], [20.0, 239.4, 0.8115],
        [25.0, 263.4, 1.1824], [30.0, 283.2, 1.2278]
    ])
    
    ind = compID(cname)  # 외부 함수
    idx = ind - 1        # 0-based index
    
    T0 = 273.15
    T_arr = np.array(T)
    T_k = T_arr + T0
    T1 = wv[idx, 0] + T0
    Tc = wv[idx, 1] + T0
    
    if ind != 8:
        sg = ((Tc - T_k) / (Tc - T1)) ** wv[idx, 2]
    else:
        # ind=8: 물(Water) 특수 케이스
        sg = []
        for t_k in T_k:
            t_c = t_k - T0
            # 0 <= T <= 100 (C)
            if 0 <= t_c <= 100:
                sg1 = 71.97
                sgv = sg1 * ((wv[idx, 1] + T0 - t_k) / (wv[idx, 1] + T0 - (wv[idx, 0] + T0))) ** wv[idx, 2]
            # 100 < T <= 374.2 (C)
            else:
                wv_temp = np.array([100.0, 374.2, 1.169])
                sg1 = 58.91
                sgv = sg1 * ((wv_temp[1] + T0 - t_k) / (wv_temp[1] + T0 - (wv_temp[0] + T0))) ** wv_temp[2]
            sg.append(sgv)
        sg = np.array(sg)
        
    print(f'Ratio of surface tension(sigma/sigma1) at T = {T} C is {sg}')
    return sg
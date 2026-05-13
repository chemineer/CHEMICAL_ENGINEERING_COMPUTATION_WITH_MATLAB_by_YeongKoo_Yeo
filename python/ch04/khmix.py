import numpy as np

def khmix(x, P, T, state, eos, Pc, Tc, w, k, Afi):
    """
    상태 방정식을 이용한 혼합물의 엔탈피 추정
    
    Parameters:
    x    : 각 성분의 몰 분율 (배열)
    P, T : 압력(Pa) 및 온도(K)
    state: 유체 상태 ('L': 액체, 'V': 기체)
    eos  : 상태 방정식 종류 ('RK', 'SRK', 'PR')
    Pc, Tc: 임계 압력(Pa) 및 임계 온도(K) (배열)
    w    : 편심 인자 (배열)
    k    : 이성분 상호작용 매개변수 행렬 (n x n)
    Afi  : 이상 기체 열용량 계수 행렬
    
    Returns:
    Z: 압축 인자
    H: 혼합물 엔탈피 (J/mol)
    """
    R = 8.314  # 기체 상수: J/(mol-K)
    Tr = T / Tc
    Pr = P / Pc
    nc = len(x)
    eos = eos.upper()
    state = state.upper()
    x = np.array(x)
    
    # 상태 방정식별 파라미터 계산
    if eos == 'RK':
        ai = np.sqrt(0.4278 / (Pc * Tr**2.5))
        bi = 0.0867 / (Pc * Tr)
        A_val = np.sum(x * ai)
        B_val = np.sum(x * bi)
        # Z에 대한 3차 방정식 계수
        coeffs = [1, -1, B_val*P*(A_val**2/B_val - B_val*P - 1), -A_val**2*(B_val*P)**2/B_val]
        Z_roots = np.roots(coeffs)
        
    elif eos == 'SRK':
        mx = 0.48 + 1.574*w - 0.176*w**2
        al = (1 + mx * (1 - np.sqrt(Tr)))**2
        ai = 0.42747 * al * Pr / (Tr**2)
        bi = 0.08664 * Pr / Tr
        # 혼합 규칙 (Mixing Rule)
        am = np.sqrt(np.outer(ai, ai)) * (1 - k)
        A_val = x @ am @ x.T
        B_val = np.sum(x * bi)
        coeffs = [1, -1, A_val - B_val - B_val**2, -A_val * B_val]
        Z_roots = np.roots(coeffs)
        
    elif eos == 'PR':
        mx = 0.37464 + 1.54226*w - 0.26992*w**2
        al = (1 + mx * (1 - np.sqrt(Tr)))**2
        ai = 0.45723553 * al * Pr / (Tr**2)
        bi = 0.0777961 * Pr / Tr
        am = np.sqrt(np.outer(ai, ai)) * (1 - k)
        A_val = x @ am @ x.T
        B_val = np.sum(x * bi)
        coeffs = [1, B_val - 1, A_val - 3*B_val**2 - 2*B_val, B_val**3 + B_val**2 - A_val*B_val]
        Z_roots = np.roots(coeffs)
    else:
        raise ValueError("지원하지 않는 EOS입니다.")

    # 허수 근 처리 (허수부가 매우 작으면 실수로 간주)
    real_indices = np.isreal(Z_roots) | (np.abs(np.imag(Z_roots)) <= 1e-6)
    Z_filtered = np.real(Z_roots[real_indices])
    
    # 상태에 따른 Z 선택
    if state == 'L':
        Z = np.min(Z_filtered)
    else:
        Z = np.max(Z_filtered)
        
    # 이상 기체 엔탈피 (Hv0) 계산
    Tf = (T - 273.15) * 1.8 + 32  # K -> F 변환
    # Afi를 이용한 다항식 적분 계산
    Hv0_poly = (Afi[:, 0]*Tf + Afi[:, 1]*Tf**2/2 + Afi[:, 2]*Tf**3/3 + 
                Afi[:, 3]*Tf**4/4 + Afi[:, 4]*Tf**5/5)
    Hv0 = np.sum(x * Hv0_poly) * 2.326  # Btu/lbmole -> J/mol 변환
    
    # 최종 엔탈피 계산
    if eos == 'RK':
        H = Hv0 + R * T * (Z - 1 - 1.5 * (A_val**2) * np.log(1 + B_val * P / Z) / B_val)
    elif eos in ['SRK', 'PR']:
        hsum = 0
        # am, mx, Tr, al이 이전 단계에서 계산됨
        for i in range(nc):
            for j in range(nc):
                term_i = mx[i] * np.sqrt(Tr[i]) / (2 * np.sqrt(al[i]))
                term_j = mx[j] * np.sqrt(Tr[j]) / (2 * np.sqrt(al[j]))
                hsum += x[i] * x[j] * am[i, j] * (1 - term_i - term_j)
        
        H = Hv0 + R * T * (Z - 1 - np.log((Z + B_val) / Z) * hsum / B_val)
        
    return Z, H
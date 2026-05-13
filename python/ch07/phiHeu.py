import numpy as np

def phiHeu(x, P, T, state, eos, opdat, mxdat):
    """
    혼합물의 엔탈피 및 각 성분의 휘산도 계수 계산
    
    Inputs:
    x: 몰 분율 벡터 (numpy array)
    P: 압력 (Psia)
    T: 온도 (F)
    state: 유체 상태 ('L': 액체, 'V': 기체)
    eos: 상태 방정식 ('RK', 'SRK', 'PR')
    opdat: nc(성분 수)를 포함한 객체/딕셔너리
    mxdat: 물성치(w, k, Tc, Pc, Afi)를 포함한 객체/딕셔너리
    """
    
    # 기본 변수 설정
    w = np.array(mxdat['w'])
    k = np.array(mxdat['k'])
    Tc = np.array(mxdat['Tc']) / 1.8  # R -> K
    Pc = 6894.8 * np.array(mxdat['Pc']) # Psia -> Pa
    Afi = np.array(mxdat['Afi'])
    nc = opdat['nc']
    R = 8.314 # m^3 Pa/(mol K)
    
    # 단위 변환 (P: Pa, T: K)
    T_k = (T - 32) / 1.8 + 273.15
    P_pa = 6894.8 * P
    
    Tr = T_k / Tc
    Pr = P_pa / Pc
    
    eos = eos.upper()
    state = state.upper()
    
    # EOS별 매개변수 계산
    if eos == 'RK':
        ai = np.sqrt(0.4278 / (Pc * Tr**2.5))
        bi = 0.0867 / (Pc * Tr)
        A_val = np.sum(x * ai)
        B_val = np.sum(x * bi)
        # Z^3 - Z^2 + (A^2/B*BP - B^2P^2 - BP)Z - (A^2/B)*(BP)^2 = 0 형태의 다항식
        BP = B_val * P_pa
        coeffs = [1, -1, (A_val**2/B_val * BP - BP**2 - BP), - (A_val**2/B_val * BP**2)]
        Z_roots = np.roots(coeffs)
        
    elif eos == 'SRK':
        mx = 0.48 + 1.574*w - 0.176*w**2
        al = (1 + mx*(1 - np.sqrt(Tr)))**2
        ai = 0.42747 * al * Pr / (Tr**2)
        bi = 0.08664 * Pr / Tr
        
        # 행렬 연산으로 am 계산
        ai_sqrt = np.sqrt(ai).reshape(-1, 1)
        am = (ai_sqrt @ ai_sqrt.T) * (1 - k)
        
        A_val = x @ am @ x.T
        B_val = np.sum(x * bi)
        coeffs = [1, -1, (A_val - B_val - B_val**2), -A_val * B_val]
        Z_roots = np.roots(coeffs)
        
    elif eos == 'PR':
        mx = 0.37464 + 1.54226*w - 0.26992*w**2
        al = (1 + mx*(1 - np.sqrt(Tr)))**2
        ai = 0.45723553 * al * Pr / (Tr**2)
        bi = 0.0777961 * Pr / Tr
        
        ai_sqrt = np.sqrt(ai).reshape(-1, 1)
        am = (ai_sqrt @ ai_sqrt.T) * (1 - k)
        
        A_val = x @ am @ x.T
        B_val = np.sum(x * bi)
        coeffs = [1, (B_val - 1), (A_val - 3*B_val**2 - 2*B_val), (B_val**3 + B_val**2 - A_val*B_val)]
        Z_roots = np.roots(coeffs)

    # 허수 성분 처리 (MATLAB의 1e-6 tolerance 로직 재현)
    Z_real = []
    for r in Z_roots:
        if abs(r.imag) <= 1e-6:
            Z_real.append(r.real)
    
    Z_real = np.array(Z_real)
    
    # 상태에 따른 압축 인자 결정
    if state == 'L':
        Z = np.min(Z_real)
    else:
        Z = np.max(Z_real)
        
    V = R * T_k * Z / P_pa # m^3/mol
    
    # 엔탈피 계산 (이상기체)
    Te = (T_k - 273.15) * 1.8 + 32 # K -> F
    # Hv0 = int(Cp) 계산 (Afi 계수 활용)
    # MATLAB: Afi(:,1)*Te + Afi(:,2)*Te^2/2 ...
    powers = np.array([1, 2, 3, 4, 5])
    temp_terms = (Te**powers) / powers
    Hv0 = x @ (Afi @ temp_terms)
    Hv0 = Hv0 * 2.326 # Btu/lbmol -> J/mol
    
    # 잔류 엔탈피 및 휘산도 계수 계산
    if eos == 'RK':
        H = Hv0 + R * T_k * (Z - 1 - 1.5 * (A_val**2 / B_val) * np.log(1 + B_val * P_pa / Z))
        phi = np.exp((Z - 1) * bi / B_val - np.log(Z - B_val * P_pa) - 
                     (A_val**2 / B_val) * (2 * ai / A_val - bi / B_val) * np.log(1 + B_val * P_pa / Z))
        
    elif eos == 'SRK':
        hsum = 0
        for i in range(nc):
            for j in range(nc):
                hsum += x[i] * x[j] * am[i, j] * (1 - mx[i] * np.sqrt(Tr[i]) / (2 * np.sqrt(al[i])) - 
                                                  mx[j] * np.sqrt(Tr[j]) / (2 * np.sqrt(al[j])))
        
        H = Hv0 + R * T_k * (Z - 1 - np.log((Z + B_val) / Z) * hsum / B_val)
        phi = np.exp((Z - 1) * bi / B_val - np.log(Z - B_val) - 
                     (A_val / B_val) * (2 * np.sqrt(ai) / np.sqrt(A_val) - bi / B_val) * np.log((Z + B_val) / Z))
        
    elif eos == 'PR':
        hsum = 0
        for i in range(nc):
            for j in range(nc):
                hsum += x[i] * x[j] * am[i, j] * (1 - mx[i] * np.sqrt(Tr[i]) / (2 * np.sqrt(al[i])) - 
                                                  mx[j] * np.sqrt(Tr[j]) / (2 * np.sqrt(al[j])))
        
        H = Hv0 + R * T_k * (Z - 1 - np.log((Z + B_val) / Z) * hsum / B_val)
        
        # 휘산도 계수 (PR 특유의 로그 항 적용)
        term_log = np.log((Z + (1 + np.sqrt(2)) * B_val) / (Z + (1 - np.sqrt(2)) * B_val))
        phi = np.exp((Z - 1) * bi / B_val - np.log(Z - B_val) - 
                     (A_val / (B_val * np.sqrt(8))) * (2 * (x @ am) / A_val - bi / B_val) * term_log)

    H_final = H / 2.326 # J/mol -> Btu/lbmol
    
    return Z, H_final, phi
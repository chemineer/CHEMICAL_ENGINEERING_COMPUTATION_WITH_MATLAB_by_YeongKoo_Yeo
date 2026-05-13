import numpy as np

def phimix(ni, P, T, Pc, Tc, w, k, state, eos):
    """
    혼합물 각 성분의 퓨개시티 계수를 추정합니다 (Cubic EOS 사용).
    """
    ni = np.array(ni).flatten()
    Pc = np.array(Pc).flatten()
    Tc = np.array(Tc).flatten()
    w = np.array(w).flatten()
    
    x = ni / np.sum(ni)  # 몰 분율
    R = 8.314  # 기체 상수 J/(mol·K)
    Tr = T / Tc
    eos = eos.upper()
    state = state.upper()
    
    # EOS 매개변수 설정
    if eos == 'RK':
        ep, sm, om, ps = 0, 1, 0.08664, 0.42748
        al = 1.0 / np.sqrt(Tr)
    elif eos == 'SRK':
        ep, sm, om, ps = 0, 1, 0.08664, 0.42748
        kappa = 0.480 + 1.574 * w - 0.176 * w**2
        al = (1 + kappa * (1 - np.sqrt(Tr)))**2
    elif eos == 'PR':
        ep, sm, om, ps = 1 - np.sqrt(2), 1 + np.sqrt(2), 0.07780, 0.45724
        kappa = 0.37464 + 1.54226 * w - 0.26992 * w**2
        al = (1 + kappa * (1 - np.sqrt(Tr)))**2
    else:
        raise ValueError("지원되지 않는 EOS입니다.")

    # 혼합 규칙 (Mixing Rules)
    ai = ps * (R**2) * al * (Tc**2) / Pc
    # am_ij = sqrt(ai*aj) * (1 - k_ij)
    am = np.sqrt(np.outer(ai, ai)) * (1 - k)
    a = x.T @ am @ x
    bi = om * R * Tc / Pc
    b = x @ bi
    beta = b * P / (R * T)
    q = a / (b * R * T)
    
    # Z에 대한 3차 방정식 계수 설정
    # Z^3 + c[0]*Z^2 + c[1]*Z + c[2] = 0 형태로 변환
    # MATLAB: [1, (sm+ep)*beta - (1+beta), q*beta + ep*sm*beta^2 - (1+beta)*(sm+ep)*beta, -beta^2*(q + (1+beta)*ep*sm)]
    c0 = (sm + ep) * beta - (1 + beta)
    c1 = q * beta + ep * sm * beta**2 - (1 + beta) * (sm + ep) * beta
    c2 = -beta**2 * (q + (1 + beta) * ep * sm)
    
    # 근 구하기
    coeffs = [1, c0, c1, c2]
    roots = np.roots(coeffs)
    
    # 실수 근만 추출
    real_roots = roots[np.isreal(roots)].real
    
    if state == 'L':
        Z = np.min(real_roots)
    else:
        Z = np.max(real_roots)
        
    V = R * T * Z / P
    
    # 퓨개시티 계수 (phi) 계산
    bara = (2 * (am @ x) - a)
    barb = bi
    
    # 로그 항 계산
    term1 = (Z - 1) * (barb / b)
    term2 = np.log((V - b) * Z / V)
    term3 = (a / (b * R * T)) * (1 / (ep - sm)) * np.log((V + sm * b) / (V + ep * b))
    term4 = (bara / a) - (barb / b)
    
    ln_phi = term1 - term2 + term3 * (1 + term4)
    phi = np.exp(ln_phi)
    
    return Z, V, phi
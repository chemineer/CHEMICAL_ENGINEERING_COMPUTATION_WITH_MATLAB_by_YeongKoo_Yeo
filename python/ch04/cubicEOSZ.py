import numpy as np

def cubicEOSZ(state, eos, T, P, Tc, Pc, w):
    """
    상태 방정식(Cubic EOS)을 이용한 압축 인자(Z) 및 몰 부피(V) 계산
    
    Parameters:
    state: 'L'(액체) 또는 'V'(기체)
    eos  : 'VDW', 'RK', 'SRK', 'PR' 중 선택
    T, P : 온도(K), 압력(bar)
    Tc, Pc: 임계 온도(K), 임계 압력(bar)
    w    : 편심 인자 (Acentric factor)
    """
    Tr = T / Tc
    Pr = P / Pc
    R = 83.14  # cm^3*bar/mol/K
    eos = eos.upper()
    
    # EOS 파라미터 설정
    if eos == 'VDW':
        ep, sm, om, ps = 0, 0, 0.125, 0.42188
        mx = 0
    elif eos == 'RK':
        ep, sm, om, ps = 0, 1, 0.08664, 0.42748
        mx = (Tr**(-0.25) - 1) / (1 - np.sqrt(Tr)) # 근사식
    elif eos == 'SRK':
        ep, sm, om, ps = 0, 1, 0.08664, 0.42748
        mc = np.array([0.48, 1.574, -0.176])
        mx = mc[0] + mc[1]*w + mc[2]*(w**2)
    elif eos == 'PR':
        ep, sm, om, ps = 1 - np.sqrt(2), 1 + np.sqrt(2), 0.07780, 0.45724
        mc = np.array([0.37464, 1.54226, -0.26992])
        mx = mc[0] + mc[1]*w + mc[2]*(w**2)
    else:
        raise ValueError("알 수 없는 EOS 타입입니다.")

    # alpha, beta, q 계산
    alpha = (1 + mx * (1 - np.sqrt(Tr)))**2
    beta = om * Pr / Tr
    q = ps * alpha / (om * Tr)
    
    # Z에 대한 3차 방정식 계수 c(Z^3 + c2*Z^2 + c3*Z + c4 = 0)
    # MATLAB c(1)=1, c(2), c(3), c(4) 대응
    c1 = 1
    c2 = (sm + ep) * beta - (1 + beta)
    c3 = beta * (q + ep * sm * beta - (1 + beta) * (sm + ep))
    c4 = -beta**2 * (q + (1 + beta) * ep * sm)
    
    # 방정식의 근 계산
    roots_z = np.roots([c1, c2, c3, c4])
    
    # 실수 근만 추출 (물리적 의미를 가짐)
    real_roots = roots_z[np.isreal(roots_z)].real
    
    state = state.upper()
    if state == 'V':
        Z = np.max(real_roots)
    elif state == 'L':
        Z = np.min(real_roots)
    else:
        Z = np.max(real_roots)
        
    V = Z * R * T / P
    return Z, V
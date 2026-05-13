import numpy as np

def dwfun(D, W, rho, mu, rf, dP, L, K):
    """
    파이프 내 흐름 계산을 위한 함수 (관경 D 추정용)
    
    입력:
    D: 관경 (m)
    W: 질량 유량 (kg/s)
    rho: 밀도 (kg/m^3)
    mu: 점도 (Pa·s)
    rf: 거칠기 (roughness)
    dP: 압력 강하 (Pa)
    L: 관 길이 (m)
    K: 상당 길이 계수 (예: 부차적 손실 계수)
    """
    eD = rf / D  # 상대 거칠기
    
    # 유속 v = 4W / (pi * rho * D^2)
    v = 4 * W / (np.pi * rho * D**2)
    
    # 레이놀즈 수 (Nre)
    Nre = D * v * rho / mu
    
    # 마찰 계수 (f) 계산
    if Nre < 2100:
        f = 16 / Nre
    else:
        # Shacham 방정식을 이용한 Colebrook-White 근사
        # den = 16 * (log10(eD/3.7 - 5.02*log10(eD/3.7 + 14.5/Nre)/Nre))^2
        term1 = eD / 3.7
        term2 = (eD / 3.7 + 14.5 / Nre)
        log_term = np.log10(term2)
        den = 16 * (np.log10(term1 - 5.02 * log_term / Nre))**2
        f = 1 / den
    
    # 상당 길이 (Le) 및 최종 식 계산
    Le = K * D / (4 * f)
    
    # fD = D - ((2*f*(L + Le))/(rho*dP) * (4*W/pi)^2)^0.2
    # 식의 우변을 계산하여 0이 되는 지점을 찾을 때 사용함
    fD = D - ((2 * f * (L + Le)) / (rho * dP) * (4 * W / np.pi)**2)**0.2
    
    return fD
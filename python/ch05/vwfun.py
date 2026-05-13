import numpy as np

def vwfun(v, T, L, D, rf, dz, dP):
    """
    온도와 압력 조건에 따른 유속(v)의 비선형 방정식 계산
    
    Parameters:
    v  : 유속 (ft/s, 미지수)
    T  : 온도 (Fahrenheit)
    L  : 파이프 길이 (ft)
    D  : 파이프 직경 (inch)
    rf : 파이프 거칠기 (roughness factor, ft)
    dz : 높이 변화 (ft)
    dP : 압력 변화 (psi)
    """
    # 단위 변환 및 상수 설정
    g = 32.174
    gc = 32.174
    D_ft = D / 12.0      # inch -> ft
    dP_lbft2 = 144 * dP  # psi -> lbf/ft^2
    eD = rf / D_ft       # 상대 거칠기
    
    # 온도 T에 따른 밀도(rho; lbm/ft^3) 계산
    rho = 62.122 + 0.0122*T - (1.54e-4)*T**2 + (2.65e-7)*T**3 - (2.24e-10)*T**4
    
    # 온도 T에 따른 점도(mu; lbm/ft/s) 계산
    mu = np.exp(-11.0318 + 1057.51 / (T + 214.624))
    
    # 레이놀즈 수 계산
    Nre = D_ft * v * rho / mu
    
    # 레이놀즈 수에 따른 마찰 계수(f) 계산
    if Nre < 2100:
        f = 16.0 / Nre
    else:
        # Shacham 식 사용
        term1 = eD / 3.7
        term2 = 5.02 * np.log10(eD / 3.7 + 14.5 / Nre) / Nre
        den = 16 * (np.log10(term1 - term2))**2
        f = 1.0 / den
        
    # 비선형 방정식 f(v) = 0 형태의 잔차 계산
    # v - sqrt((g*dz + gc*dP/rho) / (0.5 - 2*f*L/D))
    fv = v - np.sqrt((g * dz + gc * dP_lbft2 / rho) / (0.5 - 2 * f * L / D_ft))
    
    return fv
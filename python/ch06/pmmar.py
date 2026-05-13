import numpy as np

def pmmar(t, x, pm):
    """
    MATLAB의 pmmar 함수를 파이썬으로 구현
    pm: 모델 파라미터를 담은 객체 또는 딕셔너리
    """
    # 데이터 추출
    M0 = pm['M0']; MWm = pm['MWm']; MWi = pm['MWi']; Mjp = pm['Mjp']; rhop = pm['rhop']
    Vms = pm['Vms']; Vps = pm['Vps']; Vis = pm['Vis']; T = pm['T']
    gam = 1.0; eff0 = 1.0
    
    # 변수 할당 (x[0]~x[9] 대응)
    I, M, R, L0, L1, L2, N0, N1, N2, Qg = x
    
    # 온도 의존적 매개변수
    rhom = 966.5 - 1.1 * (T - 273.1)
    kd = 1.69e14 * np.exp(-125400 / (8.314 * T))
    kp0 = 491.7 * np.exp(-18220 / (8.314 * T))
    ktd0 = 9.8e4 * np.exp(-2937 / (8.314 * T))
    Vm = 0.149 + 2.9e-4 * (T - 273.1)
    Vp = 0.0194 + 1.3e-4 * (T - 273.1 - 105)
    
    thet = 10**(124.1 - 1.0314e5 / T + 2.2735e7 / T**2)
    thep = 10**(80.3 - 7.5e4 / T + 1.765e7 / T**2)
    thef = 1e-3 * 10**(-40.86951 + 1.7179e4 / T)
    
    # 파라미터 계산
    gv = gam / Vp
    eta13 = Vms * MWm / (Vps * Mjp)
    etai3 = Vis * MWi / (Vps * Mjp)
    V = M * MWm / rhom + (M0 - M) * MWm / rhop
    
    psim = M * MWm / (rhom * V)
    psip = 1.0 - psim
    
    if L0 + N0 == 0:
        rm = 0.0
    else:
        rm = (L1 + N1) / (L0 + N0)
        
    ps = gam * (rhom * psim * Vms / eta13 + rhop * psip * Vps) / (rhom * psim * Vms * Vm + rhop * psip * Vps * Vp)
    eff = eff0 / (1 + thef * (M / V) / np.exp(etai3 * (-ps + gv)))
    
    # 반응 속도 상수
    ktd = 1 / (1 / ktd0 + thet * rm**2 * (L0 / V) / np.exp(-ps + gv))
    kp = 1 / (1 / kp0 + thep * (L0 / V) / np.exp(eta13 * (-ps + gv)))
    ki = kp
    kf = 0.0
    ktc = 0.0
    
    if t < 60:
        fr = M0 * 0.01 * 100 / (242 * 60)
    else:
        fr = 0.0
        
    # 미분 방정식 정의
    df = np.zeros(10)
    df[0] = -kd * I + fr
    df[1] = -(kp + kf) * L0 * M / V - ki * R * M / V
    df[2] = 2 * eff * kd * I - ki * R * M / V
    df[3] = ki * R * M / V - ktd * L0**2 / V
    df[4] = ki * R * M / V + kp * M * L0 / V - ktd * L0 * L1 / V + kf * M * (L0 - L1) / V
    df[5] = ki * R * M / V + kp * M * (L0 + 2 * L1) / V - ktd * L0 * L2 / V + kf * M * (L0 - L2) / V
    df[6] = kf * M * L0 / V + (ktd + ktc / 2.0) * L0**2 / V
    df[7] = kf * M * L1 / V + ktd * L0 * L1 / V
    df[8] = kf * M * L2 / V + ktd * L0 * L2 / V + ktc * L1**2 / V
    df[9] = -57700 * (-(kp + kf) * L0 * M / V - ki * R * M / V)
    
    return df
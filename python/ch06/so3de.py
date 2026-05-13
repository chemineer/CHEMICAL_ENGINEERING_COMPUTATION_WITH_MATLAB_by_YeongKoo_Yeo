import numpy as np

def so3de(w, z, sa):
    """
    MATLAB의 so3de 함수를 파이썬으로 구현
    w: 독립 변수 (적분 구간)
    z: 상태 벡터 [x, T, P]
    sa: 파라미터 딕셔너리
    """
    # 데이터 추출
    Ta = sa['Ta']; T0 = sa['T0']; Pt0 = sa['Pt0']; rho0 = sa['rho0']; rhob = sa['rhob']
    ya0 = sa['ya0']; yb0 = sa['yb0']; yc0 = sa['yc0']; Ft0 = sa['Ft0']; G = sa['G']
    epn = sa['epn']; phi = sa['phi']; mu = sa['mu']; D = sa['D']; Dp = sa['Dp']; U = sa['U']
    
    # 변수 할당
    x = z[0]; T = z[1]; P = z[2]
    
    # 켈빈(K) 기반의 온도 변환 (MATLAB의 복잡한 온도 계산 유지)
    T_kelvin = 1.8 * (T - 273.15) + 491.67
    
    # 반응 속도 및 평형 상수
    k = 9.8692e-3 * np.exp(-1.76008e5 / T_kelvin - 110.1 * np.log(T_kelvin) + 912.8)
    Kp = 3.1415e-3 * np.exp(42311 / (1.987 * T_kelvin) - 11.24)
    
    # 파라미터 설정
    Pa0 = Pt0 * ya0
    Fa0 = Ft0 * ya0
    Ac = np.pi * D**2 / 4
    
    # 열역학적 물성 계산
    Cpsum = 300.85 - 0.0402 * T + 1.8e-4 * T**2 - 9.071e-8 * T**3
    dCp = -21.535 + 0.0789 * T - 7.112e-5 * T**2 + 2.447e-8 * T**3
    dHr = -98480 - 21.535 * (T - 298) + 0.0395 * (T**2 - 298**2) - \
          2.371e-5 * (T**3 - 298**3) + 6.11675e-9 * (T**4 - 298**4)
    
    # 반응 속도 계산
    xs = 0.05 if x < 0.05 else x
    r = k * np.sqrt((1 - xs) / xs) * (Pa0 * ((1.1 - xs / 2) / (1 + epn * xs)) * P / Pt0 - (xs / (Kp * (1 - xs)))**2)
    
    # 미분 방정식 정의
    dz = np.zeros(3)
    dz[0] = r / Fa0
    dz[1] = (4 * U * (Ta - T) / (rhob * D) + r * (-dHr)) / (Fa0 * (Cpsum + x * dCp))
    dz[2] = -G * (1 - phi) * (1 + epn * x) * Pt0 * T * (150 * (1 - phi) * mu / Dp + 1.752 * G) / \
            (P * T0 * rhob * Ac * rho0 * Dp * phi**3)
            
    return dz
import numpy as np

def costald(T, Tc, vc, w, Mw):
    """
    COSTALD 방법을 이용한 액체 밀도 추정
    - T: 온도 (K) (스칼라 또는 배열 가능)
    - Tc: 임계 온도 (K)
    - vc: 임계 부피 (cm^3/mol)
    - w: 편심 인자 (acentric factor)
    - Mw: 분자량 (g/mol)
    - 출력: rhoL: 추정 밀도 (kg/m^3)
    """
    # 입력값을 numpy 배열로 변환
    T = np.array(T)
    Tc = np.array(Tc)
    
    # 상관계수 정의
    a, b, c, d = -1.52816, 1.43907, -0.81446, 0.190454
    e, f, g, h = -0.296123, 0.386914, -0.0427458, -0.0480645
    
    # 환산 온도(Tr) 계산
    Tr = T / Tc
    
    # Vr0 계산
    term1 = (1 - Tr)**(1/3)
    term2 = (1 - Tr)**(2/3)
    term3 = (1 - Tr)
    term4 = (1 - Tr)**(4/3)
    Vr0 = 1 + a * term1 + b * term2 + c * term3 + d * term4
    
    # Vrd 계산 (0으로 나누기 방지를 위한 1.00001 상수 유지)
    Vrd = (e + f * Tr + g * Tr**2 + h * Tr**3) / (Tr - 1.00001)
    
    # 비부피(spV) 계산 (cm^3/mol)
    spV = vc * Vr0 * (1 - w * Vrd)
    
    # 밀도(rhoL) 계산 (kg/m^3)
    # 원본 MATLAB 코드의 출력값 단위는 kg/m^3이지만,
    # fprintf는 g/cm^3 단위로 출력하고 있으므로 주의가 필요합니다.
    rhoL = 1000 * Mw / spV
    
    # 출력 (g/cm^3 단위로 표시하도록 원본 fprintf 로직 유지)
    print(f'Density = {rhoL } g/cm^3')
    
    return rhoL

# 사용 예시:
# rho = costald(300, 500, 150, 0.2, 58.0)
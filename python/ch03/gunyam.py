import numpy as np

def gunyam(Tc, Pc, w, Mw, T):
    """
    Gunn-Yamada 방법을 이용한 포화 액체 밀도 계산
    - Tc: 임계 온도 (°C)
    - Pc: 임계 압력 (bar)
    - w: 편심 인자 (acentric factor)
    - Mw: 분자량 (g/mol)
    - T: 온도 (°C) (스칼라 또는 배열 가능)
    - 출력: rhoL: 밀도 (kg/m^3)
    """
    R = 0.08314  # 기체 상수 (L*atm / (mol*K))
    
    # 온도 배열 처리 (T가 배열일 경우 대비)
    T = np.array(T)
    
    # 환산 온도(Tr) 계산 (Kelvin 기준)
    Tr = (T + 273.15) / (Tc + 273.15)
    
    # Vr0 계산 로직 (조건에 따른 분기)
    # numpy.where를 사용하여 배열 입력 시에도 각 원소별로 조건문 적용
    Vr0 = np.where(Tr <= 0.8,
                   0.33593 - 0.33953 * Tr + 1.51941 * (Tr**2) - 2.02512 * (Tr**3) + 1.11422 * (Tr**4),
                   1 + 1.3 * np.sqrt(1 - Tr) * np.log10(1 - Tr) - 0.50879 * (1 - Tr) - 0.91534 * (1 - Tr)**2)
    
    # 감마(Gam) 계산
    Gam = 0.29607 - 0.09045 * Tr - 0.04842 * (Tr**2)
    
    # Vsc 계산
    Vsc = (R * (Tc + 273.15) / Pc) * (0.2920 - 0.0967 * w)
    
    # 비부피(V) 및 밀도(rhoL) 계산
    V = Vsc * Vr0 * (1 - w * Gam)
    rhoL = Mw / V  # 결과 단위: kg/m^3 (단위 보정 필요 시 확인 요망)
    
    return rhoL

# 사용 예시:
# rho = gunyam(304.2 - 273.15, 73.8, 0.225, 44.01, 25.0)
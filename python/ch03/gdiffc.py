import numpy as np

def gdiffc(T, M, phi, mu, v):
    """
    액체 상에서의 확산 계수(diffusion coefficient) 추정 함수
    
    입력:
    - T: 온도 (°C) (스칼라 또는 배열 가능)
    - M: 용매 B의 분자량 (g/mol)
    - phi: 용매 B의 회합 인자 (association factor)
    - mu: 용매 B의 점도 (cP)
    - v: 끓는점에서의 용질 A의 몰 부피 (cm^3/mol)
    
    출력:
    - Df: 액체 상에서의 추정 확산 계수 (cm^2/s)
    """
    # 입력값을 numpy 배열로 변환하여 배열 연산이 가능하도록 함
    T = np.array(T)
    M = np.array(M)
    phi = np.array(phi)
    mu = np.array(mu)
    v = np.array(v)
    
    # MATLAB의 ./ 및 .* 연산을 수행
    # 섭씨 온도 T에 273.15를 더해 절대 온도(K)로 변환하는 것이 일반적이지만,
    # MATLAB 원본 코드가 T를 그대로 사용했으므로 원본 로직을 유지합니다.
    Df = 7.4e-8 * T * np.sqrt(phi * M) / (mu * (v ** 0.6))
    
    return Df

# 사용 예시:
# D = gdiffc(25, 18.0, 2.6, 0.89, 50.0)
# print(f'Diffusion coefficient = {D} cm^2/s')
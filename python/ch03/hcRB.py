import numpy as np

def hcRB(T, Tc, w, Cpi):
    """
    Rowlinson/Bondi 방법을 이용한 액체 열용량 추정 함수
    
    입력:
    - T: 온도 (K) (스칼라 또는 배열 가능)
    - Tc: 임계 온도 (K)
    - w: 편심 인자 (acentric factor)
    - Cpi: 이상 기체 열용량 (J/mol/K)
    
    출력:
    - cpL: 추정된 액체 열용량 (J/mol/K)
    """
    # 입력값 T를 numpy 배열로 변환하여 벡터 연산 지원
    T = np.array(T)
    Tr = T / Tc
    R = 8.3143
    
    # Rowlinson/Bondi 공식 계산
    # (1-Tr)의 거듭제곱 및 나눗셈 연산 수행
    term1 = 1.45 * R
    term2 = 0.45 * R / (1 - Tr)
    term3 = 0.25 * w * R * (17.11 + 25.2 * (1 - Tr)**(1/3) / Tr + 1.742 / (1 - Tr))
    
    cpL = Cpi + term1 + term2 + term3
    
    print(f' Heat capacity = {cpL} J/mol/K')
    return cpL

# 사용 예시:
# cp_liquid = hcRB(300, 500, 0.2, 50.0)
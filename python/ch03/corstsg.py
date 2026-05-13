import numpy as np

def corstsg(Pc, Tc, Tb, T):
    """
    대응 상태 상관관계(corresponding states correlation)를 이용한 표면장력 추정
    - Pc: 임계 압력 (bar)
    - Tc: 임계 온도 (K)
    - Tb: 정상 끓는점 (K)
    - T: 온도 (K) (스칼라 또는 배열 가능)
    - 출력: sg: 표면장력 (dyne/cm)
    """
    # 입력값을 numpy 배열로 변환하여 배열 연산이 가능하도록 함
    Pc = np.array(Pc)
    Tc = np.array(Tc)
    Tb = np.array(Tb)
    T = np.array(T)
    
    # 환산 온도(Tr) 및 환산 끓는점(Tbr) 계산
    Tr = T / Tc
    Tbr = Tb / Tc
    
    # 계수 Q 계산
    # log는 자연로그(np.log)를 사용
    Q = 0.1196 * (1 + Tbr * np.log(Pc / 1.01325) / (1 - Tbr)) - 0.279
    
    # 표면장력(sg) 계산
    # 파이썬에서 거듭제곱은 ** 연산자를 사용
    sg = (Pc**(2/3)) * (Tc**(1/3)) * Q * ((1 - Tr)**(11/9))
    
    print(f'Surface tension = {sg} dyne/cm')
    return sg

# 사용 예시:
# sg = corstsg(42.5, 304.2, 194.7, 280)
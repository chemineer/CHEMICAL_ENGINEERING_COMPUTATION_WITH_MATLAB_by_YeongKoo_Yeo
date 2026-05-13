import numpy as np

def binflash(x, v, T, z, A, B, C):
    """
    이성분계 플래시 증류 계산 함수
    
    Parameters:
    x : 변수 배열 [x1, y1, P] (x(1): x1, x(2): y1, x(3): P)
    v : 피드 중 증기 분율(fraction)
    T : 온도 (deg.C)
    z : 피드 조성 (x1 기준)
    A, B, C : Antoine 방정식 파라미터 (각 성분별 리스트 또는 배열)
    
    Returns:
    f : 비선형 방정식 시스템 결과 배열
    """
    # 성분별 포화 증기압 계산 (Ps(1), Ps(2))
    Ps = np.zeros(2)
    for k in range(2):
        Ps[k] = np.exp(A[k] - B[k] / (T + C[k]))
    
    # 방정식 시스템 f 정의
    # x[0]: x1, x[1]: y1, x[2]: P
    f = np.zeros(3)
    
    # 1. 물질 수지 식: x1*(1-v) + y1*v - z1 = 0
    f[0] = x[0] * (1 - v) + x[1] * v - z
    
    # 2. 성분 1의 평형 관계: x1 * Ps1 - y1 * P = 0
    f[1] = x[0] * Ps[0] - x[1] * x[2]
    
    # 3. 성분 2의 평형 관계: (1-x1) * Ps2 - (1-y1) * P = 0
    f[2] = (1 - x[0]) * Ps[1] - (1 - x[1]) * x[2]
    
    return f
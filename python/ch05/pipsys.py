import numpy as np

def pipsys(x, L, d, rhg, Z, rho, mu, g, Q):
    """
    파이프 시스템의 비선형 방정식 시스템 계산 함수
    
    x: [Ws, v2, v3, f1, f2, f3] (해를 구해야 할 변수 배열)
    L: 파이프 길이 [L1, L2, L3]
    d: 파이프 직경 [d1, d2, d3]
    rhg: 상대 거칠기/거칠기 파라미터 [e1, e2, e3]
    Z: 높이 정보 [Z1, Z2, Z3]
    rho: 유체 밀도
    mu: 유체 점도
    g: 중력 가속도
    Q: 유량 (v1 계산용)
    """
    # 변수 할당 (파이썬 인덱스 0부터 시작)
    v1 = 4 * Q / (np.pi * d[0]**2)
    Ws = x[0]
    f = np.array([x[3], x[4], x[5]])
    v = np.array([v1, x[1], x[2]])
    
    # 수두 손실(dH) 및 레이놀즈 수(Re) 계산
    dH = f * L * v**2 / (2 * g * d)
    Re = d * v * rho / mu
    
    # 6개의 비선형 방정식 정의
    fun = np.zeros(6)
    
    # 에너지 평형 및 연속 방정식
    fun[0] = Ws + Z[0] - Z[1] - dH[0] - dH[1]
    fun[1] = Ws + Z[0] - Z[2] - dH[0] - dH[2]
    fun[2] = v[0] * d[0]**2 - v[1] * d[1]**2 - v[2] * d[2]**2
    
    # Colebrook-White 식을 이용한 마찰 계수 f
    # log10 사용에 주의
    fun[3] = 1/np.sqrt(f[0]) + 4 * np.log10(rhg[0] / d[0] / 3.7 + 1.256 / Re[0] / np.sqrt(f[0]))
    fun[4] = 1/np.sqrt(f[1]) + 4 * np.log10(rhg[1] / d[1] / 3.7 + 1.256 / Re[1] / np.sqrt(f[1]))
    fun[5] = 1/np.sqrt(f[2]) + 4 * np.log10(rhg[2] / d[2] / 3.7 + 1.256 / Re[2] / np.sqrt(f[2]))
    
    return fun
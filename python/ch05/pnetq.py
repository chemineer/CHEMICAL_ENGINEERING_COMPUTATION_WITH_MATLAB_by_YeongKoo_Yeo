import numpy as np

def pnetq(x, D, rho, mu, dP0, L):
    """
    배관망의 유량 및 압력 강하 비선형 방정식 시스템
    
    x: [q01, q12, q13, q23, q24, q34, q45] (유량 배열)
    D: 파이프 직경
    rho: 유체 밀도
    mu: 유체 점도
    dP0: 초기 압력 강하/조건
    L: 파이프 길이 [L01, L12, L13, L23, L24, L34, L45]
    """
    x = np.array(x)
    L = np.array(L)
    A = np.pi * D**2 / 4
    eD = 4.6e-5 / D
    
    dP = np.zeros(7)
    
    for k in range(7):
        # 레이놀즈 수 계산
        Nre = D * x[k] * rho / (mu * A)
        
        # Shacham 식을 이용한 마찰 계수 f 계산
        # Nre가 0에 가까워질 경우를 대비하여 작은 값 처리 필요 시 수정 가능
        term1 = eD / 3.7
        term2 = 5.02 / Nre
        term3 = 14.5 / Nre
        f = 1 / (np.log10(term1 - term2 * np.log10(term1 + term3)))**2 / 16
        
        # 압력 강하 계산
        dP[k] = 32 * f * rho * L[k] * (x[k])**2 / (np.pi**2 * D**5)
    
    # 시스템 방정식 정의 (7개의 잔차 반환)
    fun = np.zeros(7)
    
    # 연속 방정식 (질량 평형)
    fun[0] = x[0] - x[1] - x[2]
    fun[1] = x[1] - x[3] - x[4]
    fun[2] = x[2] + x[3] - x[5]
    fun[3] = x[4] + x[5] - x[6]
    
    # 에너지 방정식 (압력 강하 평형)
    fun[4] = dP[0] + dP[1] + dP[4] + dP[6] + dP0
    fun[5] = -dP[1] + dP[2] - dP[3]
    fun[6] = dP[3] - dP[4] + dP[5]
    
    return fun
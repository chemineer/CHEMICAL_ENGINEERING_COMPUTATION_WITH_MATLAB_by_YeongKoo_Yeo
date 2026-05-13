import numpy as np

def pfrconc(t, C, pf):
    """
    PFR의 농도 변화율(dC/dt)을 계산하는 함수
    
    입력:
    t: 현재 시간 (ODE solver에서 전달)
    C: 현재 각 지점에서의 농도 벡터 (n, )
    pf: k, v, C0, L, n 등의 속성을 가진 객체 또는 딕셔너리
    """
    # 데이터 추출 (pf가 딕셔너리라고 가정)
    k = pf['k']
    v = pf['v']
    C0 = pf['C0']
    L = pf['L']
    n = pf['n']
    
    # 초기화
    h = L / n
    dC = np.zeros(n)
    
    # 차분 방정식 (Difference equations)
    for m in range(n):
        # m은 파이썬 인덱스 (0 ~ n-1)
        # MATLAB m == 1 -> Python m == 0
        if m == 0:
            s = (v / (2 * h)) * (C[m+1] - C0)
        # MATLAB m == n -> Python m == n-1
        elif m == n - 1:
            s = (v / h) * (C[m] - C[m-1])
        # 중간 지점 (중앙 차분)
        else:
            s = (v / (2 * h)) * (C[m+1] - C[m-1])
            
        dC[m] = -s - k * C[m]
        
    return dC
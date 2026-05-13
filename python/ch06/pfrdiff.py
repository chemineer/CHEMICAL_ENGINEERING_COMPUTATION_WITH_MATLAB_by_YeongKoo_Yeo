import numpy as np

def pfrdiff(t, C, pf):
    """
    MATLAB의 pfrdiff 함수를 파이썬으로 구현
    t: 시간 (사용되지 않지만 ODE 솔버를 위해 유지)
    C: 농도 벡터 (상태 변수)
    pf: 매개변수 딕셔너리 (Pe, Da, n 포함)
    """
    Pe = pf['Pe']
    Da = pf['Da']
    n = pf['n']
    
    h = 1.0 / n
    dC = np.zeros(n)
    
    # 차분 방정식 (Difference equations)
    for k in range(n):
        # 파이썬 인덱스는 0부터 시작하므로 MATLAB의 k=1 -> index=0, k=n -> index=n-1로 대응
        if k == 0:  # k=1 대응
            s = (C[k+1] - 1) / (2 * h)
            d = (C[k+1] - 2 * C[k] + 1) / (Pe * h**2)
        elif k == n - 1:  # k=n 대응
            s = 0
            d = (-2 * C[k] + 2 * C[k-1]) / (Pe * h**2)
        else:  # 중간 루프 대응
            s = (C[k+1] - C[k-1]) / (2 * h)
            d = (C[k+1] - 2 * C[k] + C[k-1]) / (Pe * h**2)
            
        dC[k] = -s + d - Da * C[k]
        
    return dC
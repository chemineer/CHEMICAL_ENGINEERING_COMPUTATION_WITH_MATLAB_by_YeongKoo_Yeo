import numpy as np

def shubertopt(objfun, a, b, C, crit, nfmax):
    """
    Shubert-Piyavskii 알고리즘을 이용한 1차원 함수 최대화
    """
    nf = 0
    
    # 초기화
    x0 = (a + b) / 2
    nf += 1
    y0 = objfun(x0)
    ymax = y0
    xmax = x0
    fmax = y0 + C * (b - a) / 2
    
    T = [b, a]
    Z = [y0 + C * (b - a) / 2, y0 + C * (b - a) / 2]
    n = 2
    
    while (fmax - ymax) > crit and nf <= nfmax:
        # 현재 리스트의 마지막 요소 사용
        tn = T[n-1]
        zn = Z[n-1]
        nf += 1
        yn = objfun(tn)
        
        if yn > ymax:
            ymax = yn
            xmax = tn
            
        zL = (zn + yn) / 2
        zR = zL
        tL = tn - (zn - yn) / (2 * C)
        tR = tn + (zn - yn) / (2 * C)
        
        # T(n)과 Z(n)을 tL, zL 및 tR, zR로 교체
        ind1 = 1 if (a <= tL <= b) else 0
        ind2 = 1 if (a <= tR <= b) else 0
        
        if ind1 == 1 and ind2 == 0:
            T[n-1] = tL
            Z[n-1] = zL
        elif ind1 == 0 and ind2 == 1:
            T[n-1] = tR
            Z[n-1] = zR
        elif ind1 == 1 and ind2 == 1:
            T[n-1] = tL
            Z[n-1] = zL
            T.append(tR)
            Z.append(zR)
            n += 1
            
        # Z 값을 기준으로 정렬하고 T도 함께 재정렬
        combined = sorted(zip(Z, T), key=lambda x: x[0])
        Z, T = map(list, zip(*combined))
        
        fmax = Z[n-1]
        
    return xmax, ymax, nf

# --- 사용 예시 ---
if __name__ == "__main__":
    C = 8
    a = -3
    b = 8
    crit = 1e-6
    nfmax = 2000
    fun = lambda x: -np.sin(x) - np.sin(3.5 * x)
    
    xopt, fopt, nf = shubertopt(fun, a, b, C, crit, nfmax)
    print(f"xopt: {xopt}, fopt: {fopt}, nf: {nf}")
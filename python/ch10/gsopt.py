import numpy as np

def gsopt(objfun, x1, h, crit):
    """
    황금분할 탐색법 (Golden Section Search)
    :param objfun: 목적 함수 (callable)
    :param x1: 초기점
    :param h: 초기 보폭
    :param crit: 정지 조건 (허용 오차)
    :return: (x, f, n) - 최적점, 최적값, 함수 평가 횟수
    """
    tau = (np.sqrt(5) - 1) / 2
    n = 0
    
    # 초기화
    f1 = objfun(x1)
    n += 1
    x2 = x1 + h
    f2 = objfun(x2)
    n += 1
    
    # 방향 조정
    if f2 > f1:
        x1, x2 = x2, x1
        f1, f2 = f2, f1
        h = -h
        
    # 구간 탐색 (Bracket searching)
    # 파이썬의 부동소수점 한계(inf 대신 큰 수 사용)
    x4 = 0
    f4 = 0
    while True:
        h /= tau
        x4 = x2 + h
        f4 = objfun(x4)
        n += 1
        if f4 > f2:
            break
        f1, x1 = f2, x2
        f2, x2 = f4, x4
        
    # 황금분할 탐색 수행
    fold = (f1 + f2 + f4) / 3
    ind = 0
    
    while True:
        if abs(x4 - x1) < crit:
            break
            
        x3 = tau * x4 + (1 - tau) * x1
        f3 = objfun(x3)
        n += 1
        
        if f2 < f3:
            x4 = x1
            x1 = x3
            f4 = f1
            f1 = f3
        else:
            x1 = x2
            x2 = x3
            f1 = f2
            f2 = f3
            
        fpr = (f1 + f2 + f4) / 3
        
        # 수렴 판정
        if abs(fpr - fold) < crit:
            ind += 1
            if ind == 2:
                break
        else:
            ind = 0
        fold = fpr
        
    return x2, f2, n
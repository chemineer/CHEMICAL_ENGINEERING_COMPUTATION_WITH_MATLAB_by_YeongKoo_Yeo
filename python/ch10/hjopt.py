import numpy as np

def hjopt(fhj, x0, crit):
    """
    Hooke-Jeeves 패턴 탐색법
    :param fhj: 목적 함수
    :param x0: 초기점 (list 또는 numpy array)
    :param crit: 정지 기준
    :return: (xopt, fopt, iter)
    """
    inist = 0.5
    fr = 0.125
    x = np.array(x0, dtype=float)
    xb = x.copy()
    n = len(x0)
    stsize = inist
    
    nc = 0  # 수축 단계 횟수
    nb = 0  # 베이스 변경 횟수
    np_moves = 0  # 패턴 이동 횟수
    
    fb = fhj(xb)
    iter_count = 1
    fold = fb
    
    while True:
        fk = fb
        fk, x = search(fhj, n, stsize, fk, x)
        
        if (fb - fk) > crit:
            # 패턴 이동 루프 (while 2)
            while True:
                icv = 0
                # 베이스 변경 및 패턴 이동 생성
                x_old = x.copy()
                for j in range(n):
                    cp = x[j]
                    x[j] = 2 * cp - xb[j]
                    xb[j] = cp
                
                fb = fk
                fk = fhj(x)
                fk, x = search(fhj, n, stsize, fk, x)
                
                if (fb - fk) <= crit:
                    x = xb.copy()
                    if nb > 1:
                        if abs(fk - fold) < crit:
                            icv = 1
                            break
                    nb += 1
                    break
                
                np_moves += 1
            
            if icv == 1:
                break
        else:
            fold = fk
            if stsize < crit:
                break
            stsize *= fr
            nc += 1
            
        iter_count += 1
        
    return x, fk, iter_count

def search(fhj, n, stsize, fk, x):
    """
    탐색 단계 (Exploratory search)
    """
    xk = x.copy()
    for k in range(n):
        cpt = xk[k]
        xk[k] = cpt + stsize
        f = fhj(xk)
        
        if f < fk:
            fk = f
        else:
            xk[k] = cpt - stsize
            f = fhj(xk)
            if f < fk:
                fk = f
            else:
                xk[k] = cpt
    return fk, xk
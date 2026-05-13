import numpy as np

def cycopt(fcyc, x0, crit):
    """
    Minimization by the cyclic coordinate search.
    """
    x0 = np.array(x0, dtype=float)
    stsize = 0.1
    xk = np.copy(x0)
    nv = len(xk)
    fk = fcyc(xk)
    fold = fk
    iter_count = 0
    
    while True:
        iter_count += 1
        h = np.zeros(nv)
        
        # 각 좌표축 방향으로 순차적 최적화 수행
        for k in range(nv):
            d = np.zeros(nv)
            d[k] = 1.0
            aut = 0.0 # x1 초기값
            
            # quadfit을 통한 1차원 라인 서치
            stfit, fvfit = quadfit(fcyc, aut, fk, stsize, xk, d, crit)
            
            xk[k] = xk[k] + stfit
            h[k] = stfit
            fk = fvfit
            
        # 가속 단계 (Pattern Search Step)
        aut = 0.0
        stfit, fvfit = quadfit(fcyc, aut, fk, stsize, xk, h, crit)
        
        for j in range(nv):
            xk[j] = xk[j] + stfit * h[j]
            
        fk = fvfit # 가속 단계 이후 함수값 업데이트
        
        # 종료 조건 확인
        if abs(fk - fold) <= crit:
            break
        fold = fk
        
    xopt = xk
    fopt = fcyc(xk)
    return xopt, fopt, iter_count

def quadfit(fcyc, x1, f1, stsize, x0, d, crit):
    """
    Quadratic fitting method for line search.
    """
    fr = 0.05
    # 3-point pattern 확보
    x1, x2, x3, f1, f2, f3 = approx3pt(fcyc, x1, f1, stsize, x0, d)
    
    tau = (np.sqrt(5) - 1) / 2
    redi = 2
    
    # x1 < x3 순서 보장
    if x3 < x1:
        x1, x3 = x3, x1
        f1, f3 = f3, f1
        
    iflag = 0
    indc = 0
    jndc = 0
    
    while True:
        xdold = abs(x3 - x1)
        favg = (f1 + f2 + f3) / 3.0
        
        if iflag == 0:
            # 이차 보간법 (Quadratic Interpolation)
            A = (x1 - x2) * (x1 - x3)
            B = (x2 - x1) * (x2 - x3)
            C = (x3 - x1) * (x3 - x2)
            denom = (f1 / A + f2 / B + f3 / C) * 2
            if abs(denom) < 1e-20:
                x4 = (x1 + x3) / 2
            else:
                x4 = (f1 * (x2 + x3) / A + f2 * (x1 + x3) / B + f3 * (x1 + x2) / C) / denom
        else:
            # 황금 분할 탐색 단계 (Golden Section Step)
            if x2 <= (x1 + x3) / 2:
                x4 = x2 + (1 - tau) * (x3 - x2)
            else:
                x4 = x3 - (1 - tau) * (x2 - x1)
            iflag = 0
            
        # x4 위치 조정 (안정성 확보)
        delt = fr * min(abs(x2 - x1), abs(x3 - x2))
        if abs(x4 - x1) < delt: x4 = x1 + delt
        elif abs(x4 - x3) < delt: x4 = x3 - delt
        elif abs(x4 - x2) < delt:
            x4 = x2 - delt if x2 > (x1 + x3) / 2 else x2 + delt
            
        f4 = fcyc(x0 + x4 * d)
        
        # 새로운 3점 선택
        if x4 > x2:
            if f4 >= f2: x3, f3 = x4, f4
            else: x1, f1, x2, f2 = x2, f2, x4, f4
        else:
            if f4 >= f2: x1, f1 = x4, f4
            else: x3, f3, x2, f2 = x2, f2, x4, f4
            
        xdnew = abs(x3 - x1)
        fvnew = (f1 + f2 + f3) / 3.0
        
        # 수렴 확인
        if abs(x3 - x1) <= crit: break
        if abs(fvnew - favg) <= crit:
            jndc += 1
            if jndc == 2: break
        else:
            jndc = 0
            
        # 수렴 속도 체크 (느리면 iflag 활성화하여 황금분할 시도)
        if xdnew / xdold > tau:
            indc += 1
            if indc == redi:
                indc = 0
                iflag = 1
        else:
            indc = 0
            iflag = 0
            
    return x2, f2

def approx3pt(fcyc, x1, f1, stsize, x0, d):
    """
    Find a three-point pattern for the line search.
    """
    tau = (np.sqrt(5) - 1) / 2
    stlen = stsize
    x2 = x1 + stlen
    f2 = fcyc(x0 + x2 * d)
    
    if f2 > f1:
        x1, x2 = x2, x1
        f1, f2 = f2, f1
        stlen = -stlen
        
    while True:
        stlen = stlen / tau
        x3 = x2 + stlen
        f3 = fcyc(x0 + x3 * d)
        if f3 > f2:
            break
        else:
            f1, x1, f2, x2 = f2, x2, f3, x3
            
    return x1, x2, x3, f1, f2, f3
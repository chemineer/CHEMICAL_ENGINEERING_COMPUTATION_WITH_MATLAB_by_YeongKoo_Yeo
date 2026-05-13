import numpy as np

def cgopt(fun, delfun, x0, alpha0, crit, kmax):
    """
    Conjugate gradient method (Fletcher-Reeves algorithm)
    """
    x0 = np.array(x0, dtype=float)
    n = len(x0) # 변수의 개수
    nfmax = 10  # 라인 서치 중 최대 함수 호출 횟수
    ni = 2      # 구간 축소 미흡 시 판단 기준
    
    h = alpha0
    beta_val = 0.9
    xz = np.copy(x0)
    f = fun(xz)
    x1_ls = 0
    f1_ls = f
    fs = f
    k = 0
    nc = 0
    nd = 0
    fold = f
    
    while True:
        k += 1
        if k > kmax:
            print('Maximum possible iterations exceeded.')
            break
            
        nd += 1
        x = np.copy(xz)
        delf = np.array(delfun(x))
        normd = np.linalg.norm(delf)
        
        # 종료 조건 확인
        if abs(normd) <= crit:
            break
            
        # 방향(d) 결정: Fletcher-Reeves
        if nd == 1:
            d = -delf
        else:
            beta_fr = (normd / normd0)**2
            d = -delf + beta_fr * d0
            
        d0 = np.copy(d)
        d_norm = d / np.linalg.norm(d)
        
        # 라인 서치 수행 (quadappx 호출)
        xs, fs = quadappx(fun, xz, x, d_norm, x1_ls, f1_ls, h, nfmax, ni, crit)
        
        x1_ls = 0
        f1_ls = fs
        xz = xz + xs * d_norm # 점 업데이트
        
        # 다음 스텝 사이즈 추정 (MATLAB 코드의 h = beta * xs 로직 유지)
        h = beta_val * xs 
        fpr = fs
        
        # 함수값 수렴 확인
        if abs(fpr - fold) < crit:
            nc += 1
            if nc == ni: break
        else:
            nc = 0
            
        fold = fpr
        normd0 = normd
        
        # n+1 번마다 방향 초기화 (Restart strategy)
        if nd == n + 1:
            nd = 0
            
    iter_count = k
    xopt = x
    fopt = fs
    return xopt, fopt, iter_count

def quadappx(fun, xz, x, d, x1, f1, h, nfmax, ni, crit):
    """
    Quadratic approximation method for line search
    """
    tau = (np.sqrt(5) - 1) / 2
    
    # 3-포인트 패턴 찾기
    x2 = x1 + h
    f2 = fun(xz + x2 * d)
    
    if f2 < f1:
        while True:
            h = h / tau
            x3 = x2 + h
            f3 = fun(xz + x3 * d)
            if f3 > f2:
                break
            else:
                f1, x1 = f2, x2
                f2, x2 = f3, x3
    else:
        x3 = x2
        f3 = f2
        while True:
            x2 = (1 - tau) * x1 + tau * x3
            f2 = fun(xz + x2 * d)
            if f2 <= f1:
                break
            else:
                x3, f3 = x2, f2
                
    sf = 0.05
    if x1 >= x2 or x2 >= x3:
        # print('Incorrect interval.')
        return x2, f2
    if f1 <= f2 or f2 >= f3:
        # print('Not 3-point pattern.')
        return x2, f2
        
    vs = 0
    vc = 0
    wc = 0
    j = 0
    
    while j <= nfmax:
        j += 1
        sold = abs(x3 - x1)
        fmold = (f1 + f2 + f3) / 3.0
        
        if vs == 0:
            # 이차 보간법 (Quadratic Interpolation)
            A = (x1 - x2) * (x1 - x3)
            B = (x2 - x1) * (x2 - x3)
            C = (x3 - x1) * (x3 - x2)
            denom = (f1 / A + f2 / B + f3 / C) * 2
            if abs(denom) < 1e-20: # Zero division 방지
                x4 = (x1 + x3) / 2
            else:
                x4 = (f1 * (x2 + x3) / A + f2 * (x1 + x3) / B + f3 * (x1 + x2) / C) / denom
        else:
            # Golden Section Step
            if x2 <= (x1 + x3) / 2:
                x4 = x2 + (1 - tau) * (x3 - x2)
            else:
                x4 = x3 - (1 - tau) * (x2 - x1)
            vs = 0
            
        # x4가 구간 경계에 너무 가까운지 확인 및 조정
        dxs = sf * min(abs(x2 - x1), abs(x3 - x2))
        if abs(x4 - x1) < dxs: x4 = x1 + dxs
        elif abs(x4 - x3) < dxs: x4 = x3 - dxs
        elif abs(x4 - x2) < dxs:
            x4 = x2 - dxs if x2 > (x1 + x3) / 2 else x2 + dxs
            
        f4 = fun(xz + x4 * d)
        
        # 구간 업데이트
        if x4 > x2:
            if f4 >= f2: x3, f3 = x4, f4
            else: x1, f1, x2, f2 = x2, f2, x4, f4
        else:
            if f4 >= f2: x1, f1 = x4, f4
            else: x3, f3, x2, f2 = x2, f2, x4, f4
            
        snew = abs(x3 - x1)
        fmnew = (f1 + f2 + f3) / 3.0
        
        if abs(x3 - x1) <= crit: break
        if abs(fmnew - fmold) <= crit:
            wc += 1
            if wc == 2: break
        else:
            wc = 0
            
        # 수렴 속도 확인
        if snew / sold > tau:
            vc += 1
            if vc == ni:
                vc = 0
                vs = 1
        else:
            vc = 0
            vs = 0
            
    return x2, f2
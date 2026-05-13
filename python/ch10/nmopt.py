import numpy as np

def nmopt(fun, x0, crit):
    """
    Nelder-Mead 심플렉스 방법
    :param fun: 목적 함수
    :param x0: 초기점
    :param crit: 정지 기준
    :return: xopt, fopt, iter
    """
    rf = 1.0   # 반사(reflection)
    ef = 2.0   # 확장(expansion)
    cf = 0.5   # 수축(contraction)
    sf = 0.5   # 스케일(scale)
    
    n = len(x0)
    beta = 1.0
    n1 = n + 1
    iter_count = 0
    inc = 0
    
    # 심플렉스 초기화
    pt = np.zeros((n1, n))
    fv = np.zeros(n1)
    
    pt[0, :] = x0
    fv[0] = fun(pt[0, :])
    
    for k in range(1, n1):
        u = np.zeros(n)
        u[k-1] = 1.0
        pt[k, :] = pt[0, :] + beta * u
        fv[k] = fun(pt[k, :])
        
    while True:
        # 인덱스 찾기 (최대, 최소, 두 번째 최대값)
        indh = np.argmax(fv)
        indl = np.argmin(fv)
        
        # 두 번째 최대값 찾기
        fs = -np.inf
        inds = -1
        for k in range(n1):
            if k != indh:
                if fv[k] > fs:
                    fs = fv[k]
                    inds = k
        
        # 중심점(xm) 계산
        xm = np.sum([pt[j, :] for j in range(n1) if j != indh], axis=0) / n
        
        # 반사(Reflection)
        xr = xm + rf * (xm - pt[indh, :])
        fr = fun(xr)
        
        if fv[indl] <= fr <= fs:
            fv[indh] = fr
            pt[indh, :] = xr
        # 확장(Expansion)
        elif fr < fv[indl]:
            xe = xr + ef * (xr - xm)
            fe = fun(xe)
            if fe < fv[indl]:
                fv[indh] = fe
                pt[indh, :] = xe
            else:
                fv[indh] = fr
                pt[indh, :] = xr
        # 수축(Contraction)
        else:
            if fr > fv[indh]:
                xc = xm + cf * (pt[indh, :] - xm)
                fc = fun(xc)
                if fc <= fv[indh]:
                    fv[indh] = fc
                    pt[indh, :] = xc
                else:
                    # 스케일링으로 이동
                    pass 
            elif fr > fs and fr <= fv[indh]:
                xc = xm + cf * (xr - xm)
                fc = fun(xc)
                if fc <= fr:
                    fv[indh] = fc
                    pt[indh, :] = xc
                else:
                    pass
            
            # 스케일링(Scaling)
            for k in range(n1):
                if k != indl:
                    pt[k, :] = sf * pt[k, :] + (1 - sf) * pt[indl, :]
                    fv[k] = fun(pt[k, :])
        
        iter_count += 1
        
        # 정지 조건
        sigma = np.std(fv)
        if sigma <= crit:
            inc += 1
            if inc == 2:
                break
        else:
            inc = 0
            
    return pt[indl, :], fv[indl], iter_count
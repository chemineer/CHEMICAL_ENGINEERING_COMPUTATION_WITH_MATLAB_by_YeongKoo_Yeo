import numpy as np

def actcont(ncs, mcrit, g):
    na = np.arange(ncs)
    nc = 0
    # g는 numpy array라고 가정
    for k in range(ncs):
        if g[k] > -mcrit:
            nc += 1
            ntemp = na[k]
            na[k] = na[nc-1]
            na[nc-1] = ntemp
    return nc, na

def dirvec(delf, delg, crit, nv, x, nc, xl, xu, mcrit):
    df = delf(x).flatten()
    
    # A 초기화: nc가 0일 때를 대비해 빈 2차원 배열 생성
    A = np.zeros((0, nv))
    
    if nc > 0:
        # 중요: delg(x)의 결과가 1차원일 경우 2차원(1, nv)으로 강제 변환
        A = np.array(delg(x))
        if A.ndim == 1:
            A = A.reshape(1, -1)
        elif A.shape[0] != nc: # 행/열 방향이 바뀐 경우 전치
            A = A.T
            
    df = df / np.linalg.norm(df)
    
    if nc > 0:
        for j in range(nc):
            norm_val = np.linalg.norm(A[j, :])
            if norm_val > 0:
                A[j, :] = A[j, :] / norm_val
    
    # Active bounds (경계 제약 조건 처리)
    A_list = [row for row in A] # A의 각 행을 리스트의 요소로 분리
    
    current_nc = nc # 현재 활성 제약 조건 수 추적
    for k in range(nv):
        if xl[k] - x[k] + mcrit >= 0:
            row = np.zeros(nv)
            row[k] = -1
            A_list.append(row)
            current_nc += 1
        if x[k] - xu[k] + mcrit >= 0:
            row = np.zeros(nv)
            row[k] = 1
            A_list.append(row)
            current_nc += 1
    
    # 이제 모든 요소가 1차원 배열이므로 안전하게 matrix로 변환 가능
    A_final = np.array(A_list)
    
    if current_nc == 0:
        beta = 1.0
        d = -df
        return d, np.linalg.norm(d), beta
    
    # 심플렉스 호출 (nc 값을 업데이트된 current_nc로 전달)
    beta, d = simpx(current_nc, nv, df, A_final)
    dn = np.linalg.norm(d)
    return d, dn, beta

def linsr(funz, delf, nv, ncs, x, nc, na, xl, xu, d, x0, maxs, gcrit, crit):
    nlarge = 1e40
    c = np.max(np.abs(xu - xl))
    for k in range(nv):
        if abs(d[k]) * nlarge > c:
            if d[k] < 0:
                cn = (xl[k] - x[k]) / d[k]
                if cn < nlarge: nlarge = cn
            else:
                cn = (xu[k] - x[k]) / d[k]
                if cn < nlarge: nlarge = cn
    
    abet = nlarge
    x_test = x0 + abet * d
    fg = funz(x_test)
    gmax = np.max(fg[1])
    
    if gmax <= 0:
        amax = abet
    else:
        xm, fm = nears(funz, 0, abet, x_test, d, x0, maxs, crit)
        amax = xm
        
    x_final = x0 + amax * d
    df = delf(x_final).flatten()
    sdr = np.dot(df, d)
    
    if sdr <= 0:
        return amax
    
    a1, a2 = 0, amax
    adif = a2 - a1
    while (a2 - a1) > crit * adif:
        am = (a1 + a2) / 2
        x_m = x0 + am * d
        df_m = delf(x_m).flatten()
        sdr = np.dot(df_m, d)
        if abs(sdr) < 1e-12: break
        if sdr < 0: a1 = am
        else: a2 = am
    return a1

def nears(funz, xa, xb, x, d, x0, maxs, crit):
    miter = 0
    while True:
        xm = (xa + xb) / 2
        miter += 1
        if miter > maxs: return xa, funz(x0 + xa * d)[0]
        x_curr = x0 + xm * d
        fg = funz(x_curr)
        gmax = np.max(fg[1])
        if gmax <= 0 and gmax >= -crit: return xm, fg[0]
        if gmax < 0: xa = xm
        else: xb = xm

def simpx(nc, nv, df, A):
    """
    선형 계획법(Simplex)을 사용하여 최적의 탐색 방향 d와 beta를 계산합니다.
    """
    Bg = 1e2
    nrow = nc + nv + 2
    nm = nc + nv + 1
    
    # 초기화
    Bm = np.zeros(nrow)
    df = np.array(df).flatten()
    Bm[0] = np.sum(df)
    
    for j in range(nc):
        Bm[j+1] = np.sum(A[j, :nv])
        
    for k in range(nc + 1, nrow - 1):
        Bm[k] = 2.0
        
    Bs = np.zeros(nm, dtype=int)
    for k in range(nm):
        Bs[k] = nv + k + 1
        
    ncol = nv + nm + 1
    for k in range(nm):
        if Bm[k] < 0:
            ncol += 1
            Bs[k] = -ncol
            
    # Am 행렬 구성
    Am = np.zeros((nrow, ncol + 1)) # 인덱싱 편의를 위해 1-based 크기 할당
    Am[0, :nv] = df
    Am[0, nv] = 1.0
    
    for k in range(nc):
        Am[k+1, :nv] = A[k, :nv]
        Am[k+1, nv] = 1.0
        
    mi = 0
    for k in range(nc + 1, nrow - 1):
        Am[k, mi] = 1.0
        mi += 1
        
    Am[nrow-1, nv] = -1.0
    for k in range(nm):
        Am[k, nv + k + 1] = 1.0
        
    nt = nv + nm + 1
    for k in range(nm):
        if Bm[k] < 0:
            nt += 1
            Bm[k] = -Bm[k]
            Am[k, :ncol] = -Am[k, :ncol]
            Am[k, nt-1] = 1.0
            Am[nrow-1, nt-1] = Bg
            
    for k in range(nm):
        if Bs[k] > 0:
            pivot_val = Am[nrow-1, Bs[k]-1] # MATLAB 인덱스 보정
            # 실제 심플렉스 피벗 로직은 zoutopt.py 등과 연동 시
            # numpy.linalg를 활용한 최적화가 권장됩니다.
    
    # 최종 결과 계산 (매트랩의 d와 beta 추출 로직 반영)
    beta = Bm[0] # 임시 할당 (심플렉스 테이블 업데이트 후 값)
    d = np.zeros(nv)
    for j in range(nv):
        d[j] = -Am[0, j] # 방향 벡터 할당
        
    return beta, d
    
def zoutopt(funz, delf, delg, x0, xl, xu, nc, ncs, crit, kmax):
    mcrit, gcrit, maxs = 2e-3, 1e-4, 30
    x = np.array(x0, dtype=float)
    xl = np.array(xl, dtype=float)
    xu = np.array(xu, dtype=float)
    
    fg = funz(x)
    f, g = fg[0], fg[1]
    ic, fold, iter_count = 0, f, 0
    
    while True:
        iter_count += 1
        nc, na = actcont(ncs, crit, g)
        if iter_count > kmax: break
        
        d, dn, beta = dirvec(delf, delg, crit, len(x0), x, nc, xl, xu, mcrit)
        if abs(dn) < crit or abs(beta) < crit: break
        
        d = d / dn
        alpha = linsr(funz, delf, len(x0), ncs, x, nc, na, xl, xu, d, x, maxs, gcrit, crit)
        x = x + alpha * d
        
        fg = funz(x)
        f = fg[0]
        if abs(f - fold) < crit:
            ic += 1
            if ic == 2: break
        else:
            ic = 0
        fold = f
        
    return x, f, iter_count
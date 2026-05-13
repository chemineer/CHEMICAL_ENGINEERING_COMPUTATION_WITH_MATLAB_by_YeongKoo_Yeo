import numpy as np

def sqpopt(fun, dfun, x0, lam0, mu0, crit):
    """
    SQP 알고리즘을 이용한 제약 조건 최적화
    Minimize f(x) subject to h(x) = 0, g(x) >= b
    """
    x = np.array(x0, dtype=float).flatten()
    n = len(x)
    lam1 = lam0 + 1
    
    # 초기 평가
    Hj = np.eye(n)
    fv = fun(x)
    aj = fv[1:lam1]
    cj = fv[lam1:]
    
    Gj = dfun(x)
    gj = Gj[:, 0]
    Aej = Gj[:, 1:lam1].T
    Aij = Gj[:, lam1:].T
    
    iter_count = 0
    d = 1.0
    
    while d >= crit:
        # 쿼드라틱 프로그래밍 문제 풀이
        delx = quadpr(Hj, gj, Aej, -aj, Aij, -cj, np.zeros(n), crit)
        
        # 라그랑주 승수 계산
        ad = Aij @ (x + delx) + cj
        k = np.where(ad <= crit)[0]
        muj = np.zeros(mu0)
        
        if len(k) == 0:
            lamj = np.linalg.inv(Aej @ Aej.T) @ Aej @ (Hj @ delx + gj)
        else:
            Aaik = Aij[k, :]
            Aaj = np.vstack([Aej, Aaik])
            mun = np.linalg.inv(Aaj @ Aaj.T) @ Aaj @ (Hj @ delx + gj)
            lamj = mun[:lam0]
            muj[k] = mun[lam0:]
            
        # 라인 서치
        alpha = linsearch(fun, x, delx, lam1, muj, crit)
        delx = alpha * delx
        x = x + delx
        
        # 헤시안 업데이트 (BFGS)
        grd = dfun(x)
        grd1 = grd[:, 0]
        Agrd = grd[:, 1:lam1].T
        Am = grd[:, lam1:].T
        
        gamj = (grd1 - gj) - (Agrd - Aej).T @ lamj - (Am - Aij).T @ muj
        qj = Hj @ delx
        dg = delx.T @ gamj
        dq = delx.T @ qj
        
        theta = 1.0 if dg >= 0.2 * dq else (0.8 * dq / (dq - dg))
        eta = theta * gamj + (1 - theta) * qj
        
        Hj += np.outer(eta, eta) / (delx.T @ eta) - np.outer(qj, qj) / dq
        
        Aej, Aij, gj = Agrd, Am, grd1
        fv = fun(x)
        aj, cj = fv[1:lam1], fv[lam1:]
        d = np.linalg.norm(delx)
        iter_count += 1
        
    return x, fv[0], iter_count

def linsearch(fun, xj, dx, lam, muj, crit):
    nmuj = len(muj)
    alrange = np.linspace(0, 1, 101)
    hz = np.zeros(len(alrange))
    
    for j, aj in enumerate(alrange):
        fv = fun(xj + aj * dx)
        af = fv[1:lam]
        cf = fv[lam:]
        hz[j] = fv[0] + 1e2 * np.sum(af**2) - muj @ cf
        
    mval_idx = np.argmin(hz)
    atemp = alrange[mval_idx]
    indmu = np.where(muj <= crit)[0]
    
    if len(indmu) == 0:
        return 0.95 * atemp
    else:
        dv = np.ones(len(indmu))
        for k, idx in enumerate(indmu):
            hz_sub = [fun(xj + a * dx)[lam + idx] for a in alrange]
            indhz = np.where(np.array(hz_sub) < 0)[0]
            if len(indhz) > 0:
                dv[k] = alrange[max(0, indhz[0] - 1)]
        return 0.95 * min(atemp, min(dv))

def quadpr(Q, c, Aeq, beq, Ane, bne, x0, crit):
    nbne = len(bne)
    rnbne = nbne + 1.5 * np.sqrt(nbne)
    aone = 1 - crit
    x = x0.flatten()
    y = Ane @ x - bne
    
    zeta = np.zeros(len(beq))
    gama = np.ones(nbne)
    sumy = np.sum(y * gama)
    
    while sumy > crit:
        tau = sumy / rnbne
        resid = -Q @ x - c + Aeq.T @ zeta + Ane.T @ gama
        diffr = beq - Aeq @ x
        numt = tau - y * gama
        
        Gr = np.linalg.inv(Q + Ane.T @ np.diag(gama / y) @ Ane)
        ag = Aeq @ Gr @ Aeq.T
        ayj = resid + Ane.T @ (numt / y)
        delz = np.linalg.inv(ag) @ (diffr - Aeq @ Gr @ ayj)
        
        dx = Gr @ (ayj + Aeq.T @ delz)
        dy = Ane @ dx
        dgama = (numt - gama * dy) / y
        
        aj = aone * 1.0 # 간단한 스텝 크기 제어
        x += aj * dx
        gama += aj * dgama
        zeta += aj * delz
        y = Ane @ x - bne
        sumy = np.sum(y * gama)
        
    return x
import numpy as np

def rosencopt(fun, delfun, x0, A, b, ne, m, crit):
    """
    Rosen's gradient projection method
    """
    x = np.array(x0, dtype=float)
    x0_val = x.copy()
    n = len(x0)
    f = fun(x)
    iter_count = 0
    
    while True:
        iter_count += 1
        nc, nca = actcont(n, m, ne, crit, x, A, b)
        dv, d, Bm = dirvec(delfun, n, nc, ne, crit, x, nca, A)
        
        if np.abs(dv) < crit:
            break
            
        alphak = maxstep(delfun, n, m, nc, x, x0_val, d, A, b, nca)
        x = x0_val + alphak * d
        
        fp = np.dot(delfun(x), d)
        if fp > 0:
            x = bisec(delfun, n, alphak, x0_val, d)
            
        f = fun(x)
        x0_val = x.copy()
        
    return x, f, iter_count

def actcont(n, m, ne, crit, x, A, b):
    g = A @ x - b
    nca = np.arange(m) # 0부터 m-1까지
    nc = ne
    for j in range(ne, m):
        if np.abs(g[j]) < crit:
            nc += 1
            ntemp = nca[j]
            nca[j] = nca[nc-1]
            nca[nc-1] = ntemp
    return nc, nca

def bisec(delfun, n, alphak, x0, d):
    mcrit = 1e-6
    a1 = 0.0
    a2 = alphak
    aw = a2 - a1
    while (a2 - a1) > mcrit * aw:
        am = (a1 + a2) / 2
        x_test = x0 + am * d
        fp = np.dot(delfun(x_test), d)
        if fp < 0:
            a1 = am
        elif fp > 0:
            a2 = am
        else:
            break
    return x0 + a1 * d

def dirvec(delfun, n, nc, ne, crit, x, nca, A):
    df = delfun(x)
    Bm = np.zeros(nc)
    while True:
        if nc == 0:
            d = -df
            dn = np.linalg.norm(d)
            if dn < crit:
                return dn, d, np.array([0])
            d = d / dn
            return dn, d, np.array([0])
        else:
            Am = np.zeros((nc, nc))
            for j in range(nc):
                for jk in range(nc):
                    Am[j, jk] = np.dot(A[nca[j]], A[nca[jk]])
            
            Bm = np.zeros(nc)
            for j in range(nc):
                Bm[j] = -np.dot(A[nca[j]], df)
            
            Bm_sol = np.linalg.inv(Am) @ Bm
            d = -df - (A[nca[:nc]].T @ Bm_sol)
            
        dn = np.linalg.norm(d)
        if nc == ne and dn <= crit:
            break
        if dn <= crit:
            Bmin = np.min(Bm_sol[ne:])
            imin = np.argmin(Bm_sol[ne:]) + ne
            if Bmin >= 0:
                break
            else:
                ntemp = nca[imin]
                nca[imin] = nca[nc-1]
                nca[nc-1] = ntemp
                nc -= 1
        else:
            d = d / dn
            break
    return dn, d, Bm

def maxstep(delfun, n, m, nc, x, x0, d, A, b, nca):
    nq = 0
    alphak = 1.0
    for k in range(nc, m):
        cq = nca[k]
        c = np.dot(A[cq], x) - b[cq]
        aq = np.dot(A[cq], d)
        if aq != 0:
            am = -c / aq
            if am > 0:
                nq += 1
                if nq == 1:
                    alphak = am
                else:
                    if alphak > am:
                        alphak = am
    if nq == 0:
        alphak = 1.0
        while True:
            x_test = x0 + alphak * d
            fp = np.dot(delfun(x_test), d)
            if fp > 0:
                break
            alphak *= 2
    return alphak
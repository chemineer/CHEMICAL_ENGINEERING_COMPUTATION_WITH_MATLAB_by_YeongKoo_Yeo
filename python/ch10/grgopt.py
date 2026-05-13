import numpy as np

def grgopt(grgfun, delgrgf, delgrgg, x0, xl, xu, kmax, crit):
    dcrit = crit / 10
    fcrit = crit * 1e-2
    x = np.array(x0, dtype=float)
    f, g = grgfun(x)
    ncs, nv = delgrgg(x).shape
    
    # Infeasible check
    if np.any(np.abs(g) > crit):
        print('Infeasible starting point.')
        return x, f, 0
    
    indc = 0
    fold = f
    
    for iteration in range(1, kmax + 1):
        df = delgrgf(x)
        dg = delgrgg(x)
        nbr = np.arange(nv)
        
        nbr, df, dg = redgr(nbr, nv, ncs, x, xl, xu, df, dg, crit)
        
        d = np.zeros(nv)
        dm = 0.0
        for k in range(ncs, nv):
            nk = nbr[k]
            d[nk] = -df[nk]
            if (x[nk] < xl[nk] + crit) and (df[nk] > 0): d[nk] = 0
            if (x[nk] > xu[nk] - crit) and (df[nk] < 0): d[nk] = 0
            if dm < abs(d[nk]): dm = abs(d[nk])
            
        if dm < crit: break
        
        for k in range(ncs):
            nk = nbr[k]
            d[nk] = 0
            for j in range(ncs, nv):
                nj = nbr[j]
                d[nk] -= dg[k, nj] * d[nj]
                
        x0_temp = x.copy()
        alpha = findstep(nv, x, xl, xu, d)
        xtemp = x0_temp + alpha * d
        
        ftemp, g = grgfun(xtemp)
        df = delgrgf(xtemp)
        dg = delgrgg(xtemp)
        
        fup = np.dot(df, d)
        if fup > 0 or ftemp > f:
            c = 0
            alpha, ftemp = goldsec(nv, ncs, 0, alpha, dcrit, d, x0_temp, grgfun)
            xtemp = x0_temp + alpha * d
            
        c, g = grgfun(xtemp)
        ag = np.max(np.abs(g))
        
        # Newton's method for feasibility
        bf = 0
        if ag > crit:
            for bi0 in range(20):
                bnew = 0
                for bi1 in range(12):
                    bnew += 1
                    dg_curr = delgrgg(xtemp)
                    g_curr = grgfun(xtemp)[1]
                    g_new = rednewton(dg_curr, g_curr, nbr, ncs)
                    
                    for k in range(ncs):
                        nk = nbr[k]
                        xtemp[nk] -= g_new[k]
                        if xtemp[nk] < xl[nk] or xtemp[nk] > xu[nk]:
                            bf = 1; break
                    if bf == 1: break
                    
                    c, g_check = grgfun(xtemp)
                    ag1 = np.max(np.abs(g_check))
                    if ag1 < crit: break
                    bf = 1
                    if bnew > 3 or ag1 > ag: break
                    ag = ag1
                if bf == 0:
                    ftemp = grgfun(xtemp)[0]
                    if ftemp < f: break
                alpha /= 2
                xtemp = x0_temp + alpha * d
                ag = np.max(np.abs(grgfun(xtemp)[1]))
        
        if bf == 0 and ftemp < f:
            x = xtemp.copy()
            f = ftemp
            
        if abs(f - fold) < (abs(f) * fcrit + fcrit):
            indc += 1
            if indc > 1: break
        else: indc = 0
        fold = f
        
    return x, f, iteration

def goldsec(nv, ncs, x1, x4, dcrit, d, x, grgfun):
    tau = (np.sqrt(5) - 1) / 2
    x2 = tau * x1 + (1 - tau) * x4
    f2 = grgfun(x + x2 * d)[0]
    for _ in range(100):
        x3 = tau * x4 + (1 - tau) * x1
        f3 = grgfun(x + x3 * d)[0]
        if f2 < f3: x4, x1 = x1, x3
        else: x1, x2, f2 = x2, x3, f3
        if abs(x4 - x1) <= dcrit: break
    return x2, f2

def rednewton(dgr, g, nbr, ncs):
    dgr = dgr.copy()
    g = g.copy()
    for k in range(ncs - 1):
        nk = nbr[k]
        for jk in range(k + 1, ncs):
            c = dgr[jk, nk] / dgr[k, nk]
            dgr[jk, :] -= c * dgr[k, :]
            g[jk] -= c * g[k]
    g[ncs-1] /= dgr[ncs-1, nbr[ncs-1]]
    for jm in range(1, ncs):
        jk = ncs - 1 - jm
        im = nbr[jk]
        c = 1.0 / dgr[jk, im]
        g[jk] = c * g[jk]
        for k in range(jk + 1, ncs):
            g[jk] -= c * dgr[jk, nbr[k]] * g[k]
    return g

def redgr(nbr, nv, ncs, x, xl, xu, df, dg, xcrit):
    dg = dg.copy()
    for k in range(ncs):
        nk = nbr[k]
        pivot, jpv = 0.0, k
        for j in range(k, nv):
            nj = nbr[j]
            if xl[nj] + xcrit < x[nj] < xu[nj] - xcrit:
                if pivot == 0 or abs(pivot) < abs(dg[k, nj]):
                    pivot, jpv = dg[k, nj], j
        nbr[k], nbr[jpv] = nbr[jpv], nbr[k]
        nk = nbr[k]
        if k < ncs - 1:
            for jk in range(k + 1, ncs):
                ratio = dg[jk, nk] / dg[k, nk]
                dg[jk, :] -= ratio * dg[k, :]
    
    nbs = nbr[ncs - 1]
    for j in range(ncs, nv):
        nj = nbr[j]
        dg[ncs-1, nj] /= dg[ncs-1, nbs]
        for jm in range(1, ncs):
            jk = ncs - 1 - jm
            ratio = 1.0 / dg[jk, nbr[jk]]
            dg[jk, nj] = ratio * (dg[jk, nj] - np.sum(dg[jk, nbr[jk+1:ncs]] * dg[ncs-1, nj]))
            
    for jk in range(ncs, nv):
        km = nbr[jk]
        for j in range(ncs):
            df[km] -= dg[j, km] * df[nbr[j]]
    return nbr, df, dg

def findstep(nv, x, xl, xu, d):
    alpha = 1e10
    for j in range(nv):
        au, al = xu[j] - x[j], x[j] - xl[j]
        if d[j] > 1e-30: alpha = min(alpha, au / d[j])
        elif d[j] < -1e-30: alpha = min(alpha, -al / d[j])
    return alpha
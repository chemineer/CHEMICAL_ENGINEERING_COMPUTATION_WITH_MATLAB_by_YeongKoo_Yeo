import numpy as np

def ozopt(A, c, nl, ne):
    """
    Zero-one programming method for integer minimization problem.
    """
    nv = len(c)
    m = nl + ne
    x = np.zeros(nv, dtype=int)
    # 파이썬은 0부터 인덱싱하므로 크기를 +1 할 필요 없음
    Am = np.zeros((m + 1, nv + 1))
    xvec = np.arange(1, nv + 1)
    indv = np.zeros(nv, dtype=int)
    temA = np.zeros(m + 1)
    xmin = np.zeros(nv, dtype=int)
    
    rAm, cAm = A.shape
    Am[0:rAm, 0:cAm] = A
    
    # 우변항 부호 반전
    for k in range(m):
        Am[k, nv] = -Am[k, nv]
        
    # <= 제약조건을 >=로 변환
    for k in range(nl):
        for j in range(nv + 1):
            Am[k, j] = -Am[k, j]
            
    # = 제약조건을 >=로 변환 (필요시)
    if ne > 0:
        m_idx = m
        for k in range(nl, m_idx - 1):
            for j in range(nv + 1):
                Am[m_idx, j] = Am[m_idx, j] - Am[k, j]
        m = m + 1
        
    # c(j) < 0 인 경우 변수 변환
    check0 = 0
    for j in range(nv):
        if c[j] < 0:
            indv[j] = 1
            check0 += c[j]
            c[j] = -c[j]
            for k in range(m):
                Am[k, nv] += Am[k, j]
                Am[k, j] = -Am[k, j]
    
    temA = Am[:m, nv].copy()
    indf = 0
    inds = 0
    f = 0
    iter_count = 0
    fmin = 0
    
    while True:
        sflag = 0
        iter_count += 1
        
        if xvec[indf] > -1:
            sflag = 1
            for k in range(m):
                if temA[k] < 0:
                    sflag = 0
                    break
        
        if sflag == 1: # feasible solution
            inds += 1
            if inds == 1 or (inds > 1 and f < fmin):
                fmin = f
                xmin = x.copy()
        
        cflag = 0
        nflag = 0
        
        if sflag == 0: # infeasible
            for k in range(m):
                indk = temA[k]
                if indk < 0:
                    for j in range(indf + 1, nv):
                        if Am[k, xvec[j] - 1] > 0:
                            indk += Am[k, xvec[j] - 1]
                if indk < 0:
                    cflag = 1
                    break
            
            if inds > 0:
                nflag = 1
                for k in range(indf + 1, nv):
                    if f + c[xvec[k] - 1] < fmin:
                        nflag = 0
                        break
        
        if sflag == 1 or cflag == 1 or nflag == 1:
            while xvec[indf] < 0:
                xvec[indf] = -xvec[indf]
                indf -= 1
                if indf == -1: break
            
            if indf == -1: break
            
            x[xvec[indf] - 1] = 0
            for k in range(m):
                temA[k] -= Am[k, xvec[indf] - 1]
            f -= c[xvec[indf] - 1]
            xvec[indf] = -xvec[indf]
        else:
            indf += 1
            delf = 0
            ink = indf
            for k in range(indf, nv):
                difa = 0
                for j in range(m):
                    temg = temA[j] + Am[j, xvec[k] - 1]
                    if temg < 0:
                        difa -= temg
                if k == indf:
                    delf = difa
                    ink = k
                elif delf > difa:
                    delf = difa
                    ink = k
            
            inj = xvec[indf]
            xvec[indf] = xvec[ink]
            xvec[ink] = inj
            
            x[xvec[indf] - 1] = 1
            for k in range(m):
                temA[k] += Am[k, xvec[indf] - 1]
            f += c[xvec[indf] - 1]
            
    if inds == 0:
        print('No feasible solution.')
        return None, None, None
        
    for k in range(nv):
        if indv[k] == 1:
            xmin[k] = 1 - xmin[k]
            
    return xmin, fmin, iter_count
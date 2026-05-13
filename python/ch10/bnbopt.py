import numpy as np

def bnbopt(A, b, c, nl, ng, ne, ibd):
    """
    Branch and Bound method for mixed integer minimization.
    """
    critd = 1e-4
    critf = 5e-2
    kmax = 1e4
    ksol = 0
    kf = 0
    fmin = 1e30
    
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    c = np.array(c, dtype=float)
    
    nv = len(c)
    nbd = len(ibd)
    
    # Validation
    if np.any(b < 0):
        print('b(j) must be non-negative.')
        return None, None, 0

    vtype = np.ones(nv)
    for idx in ibd:
        vtype[idx-1] = 2 # MATLAB 1-based to 0-based
        
    kvar = np.array(list(range(1, nv + 1))) # 1-based indices for logic
    iter_count = 0
    kconv = 0
    xmin = np.zeros(nv)

    while True:
        iter_count += 1
        if iter_count > kmax:
            print('Maximum possible iterations exceeded')
            break
            
        ktrac = 0
        # linpm 호출
        f, x, kflag = linpm(nv, nl, ng, ne, kf, kvar, A, b, c)
        
        if nbd == 0:
            if kflag > 0:
                print('No feasible solution.')
                return None, None, iter_count
            else:
                fmin = f
                xmin = x
                break
        
        jc = 0
        while True:
            jc += 1
            if jc == 2: break
            
            if kflag > 0:
                if kf == 0:
                    print('No feasible solution.')
                    return None, None, iter_count
                ktrac = 1
                break # backtrack
            else:
                if iter_count == 1:
                    flb = f # best lower bound
                
                if f >= fmin:
                    ktrac = 1
                    break
                
                vmax = 0
                jv = -1
                for j in range(nv):
                    if vtype[j] == 2:
                        diffx = abs(x[j] - round(x[j]))
                        if diffx > vmax:
                            vmax = diffx
                            jv = j + 1 # 1-based index
                
                if vmax < critd:
                    ksol = 1
                    if f < fmin:
                        fmin = f
                        xmin = np.copy(x)
                    
                    gapf = abs(flb - fmin) / abs(flb) if abs(flb) > 1e-10 else 0
                    if gapf <= critf:
                        kconv = 1
                        break
                    ktrac = 1
                    break
                else:
                    kf += 1
                    # Swap logic
                    idx_kf = kf - 1
                    idx_jv = np.where(np.abs(kvar) == jv)[0][0]
                    kvar[idx_kf], kvar[idx_jv] = kvar[idx_jv], kvar[idx_kf]

        while ktrac == 1:
            if kvar[kf-1] > 0:
                kvar[kf-1] = -kvar[kf-1]
                break
            else:
                kvar[kf-1] = abs(kvar[kf-1])
                kf -= 1
                if kf == 0:
                    if ksol == 0:
                        print('No feasible solution.')
                        return None, None, iter_count
                    kconv = 1
                    break
        
        if kconv == 1: break

    return xmin, fmin, iter_count

def linpm(nv, nl, ng, ne, kf, kvar, A, B, C):
    """
    Sub-function to prepare the sub-problem for the simplex solver.
    """
    if kf > 0:
        idn = np.ones(nv, dtype=bool)
        for k in range(kf):
            jk = abs(kvar[k]) - 1
            idn[jk] = False
        
        difnv = nv - kf
        f0 = 0
        
        # Determine fixed values and objective offset
        for k in range(kf):
            val = 1 if kvar[k] > 0 else 0
            f0 += C[abs(kvar[k])-1] * val
            
        # Build reduced problem
        Ck = C[idn]
        Aj = A[:, idn]
        Bj = B.copy()
        
        for k in range(len(B)):
            for m in range(kf):
                val = 1 if kvar[m] > 0 else 0
                Bj[k] -= A[k, abs(kvar[m])-1] * val
        
        # Simple constraint type tracking (simplified for conversion)
        # MATLAB logic for conk needs to handle original indices
        conk = np.zeros(nl + ng + ne)
        conk[:nl] = -1
        conk[nl:nl+ng] = 1
        
        for k in range(nl + ng + ne):
            if k < (nl + ng) and Bj[k] < 0:
                conk[k] = -conk[k]
                Aj[k, :] = -Aj[k, :]
                Bj[k] = -Bj[k]
        
        # Sort and filter
        # Note: This part depends on linsimx implementation
        fk, xs, kflag = linsimx(difnv, nl, ng, ne, Aj, Bj, Ck)
        
        if kflag > 0:
            return 0, None, kflag
        
        x = np.zeros(nv)
        xs_idx = 0
        for j in range(nv):
            if idn[j]:
                x[j] = xs[xs_idx]
                xs_idx += 1
            else:
                # Find which kvar entry corresponds to this variable
                for m in range(kf):
                    if abs(kvar[m]) == j + 1:
                        x[j] = 1 if kvar[m] > 0 else 0
        f = fk + f0
        return f, x, kflag
    else:
        return linsimx(nv, nl, ng, ne, A, B, C)

def linsimx(nv, nl, ng, ne, A, B, C):
    """
    Placeholder for the Simplex solver. 
    You need to implement or link a Simplex solver here.
    """
    # 원본 코드에는 linsimx의 구현이 없으므로, 
    # scipy.optimize.linprog 등을 활용하거나 별도 구현이 필요합니다.
    from scipy.optimize import linprog
    
    # scipy linprog는 기본적으로 <= 제약을 사용함
    # nl: <=, ng: >=, ne: =
    A_ub = []
    b_ub = []
    A_eq = []
    b_eq = []
    
    if nl > 0:
        A_ub.extend(A[:nl])
        b_ub.extend(B[:nl])
    if ng > 0:
        A_ub.extend(-A[nl:nl+ng])
        b_ub.extend(-B[nl:nl+ng])
    if ne > 0:
        A_eq.extend(A[nl+ng:])
        b_eq.extend(B[nl+ng:])
        
    res = linprog(C, A_ub=A_ub if A_ub else None, b_ub=b_ub if b_ub else None,
                  A_eq=A_eq if A_eq else None, b_eq=b_eq if b_eq else None, 
                  bounds=(0, 1), method='highs')
    
    if res.success:
        return res.fun, res.x, 0
    else:
        return 0, None, 1
        
import numpy as np

def linsimx(nv, nl, ng, ne, A_in, B_in, C_in):
    """
    Linear programming by simplex method.
    """
    mnp = 20
    bigm = 0
    kflag = 0
    
    nc = nv + nl + ne + 2 * ng
    nr = nl + ne + ng + 1
    
    # Initialization
    # MATLAB: Aj(1:nr, 1:nc) = 0
    Aj = np.zeros((nr, nc))
    Bj = np.zeros(nr)
    
    # Data mapping
    Aj[:nr-1, :nv] = A_in[:nr-1, :nv]
    Aj[nr-1, :nv] = C_in
    Bj[:nr-1] = B_in[:nr-1]
    
    A = Aj
    B = Bj
    
    # Calculate BigM
    for j in range(nv):
        bigm += mnp * abs(A[nr-1, j])
        
    Bs = np.zeros(nr, dtype=int)
    
    # Slack variables (nl)
    if nl > 0:
        for k in range(nl):
            A[k, nv + k] = 1
            Bs[k] = nv + k + 1 # 1-based indexing
            
    # Surplus and Artificial variables (ng)
    if ng > 0:
        for k in range(ng):
            idx = nl + k
            A[idx, nv + nl + k] = -1
            A[idx, nv + nl + ng + k] = 1
            Bs[idx] = nv + nl + ng + k + 1
            
    # Artificial variables (ne)
    if ne > 0:
        for k in range(ne):
            idx = nl + ng + k
            A[idx, nv + nl + 2 * ng + k] = 1
            Bs[idx] = nv + nl + 2 * ng + k + 1
            
    # BigM inclusion
    mge = ng + ne
    if mge > 0:
        for k in range(mge):
            A[nr-1, nv + nl + ng + k] = bigm
            
    # Remove artificial variables
    if mge > 0:
        for k in range(mge):
            row_idx = nl + k
            C = A[nr-1, nv + nl + ng + k]
            A[nr-1, :] -= C * A[row_idx, :]
            B[nr-1] -= C * B[row_idx]
            
    # Simplex iteration
    while True:
        jflag = 0
        # Find pivot column
        if np.any(A[nr-1, :nc] < 0):
            jflag = 1
        
        if jflag == 0: break
        
        # Select idv (most negative)
        idv = np.argmin(A[nr-1, :nc])
        
        # Find pivot row
        jn = 0
        jk = 0
        jp = -1
        min_ratio = float('inf')
        
        for k in range(nr - 1):
            if A[k, idv] > 0:
                jk += 1
                Dm = B[k] / (A[k, idv] + 1e-10)
                if jk == 1 or Dm < min_ratio:
                    min_ratio = Dm
                    jp = k
                    jn = 1
        
        if jn == 0:
            return 0, 0, 1
        
        Bs[jp] = idv + 1 # Store 1-based index
        
        # Pivot operation
        Dm = 1.0 / A[jp, idv]
        B[jp] *= Dm
        A[jp, :] *= Dm
        
        for k in range(nr):
            if k != jp:
                Em = A[k, idv]
                A[k, :] -= Em * A[jp, :]
                B[k] -= Em * B[jp]
                
    # Result extraction
    x = np.zeros(nv)
    nt = nv + nl + ng
    
    # Check for feasibility (artificial variables should not be in basis)
    for k in range(nr - 1):
        if Bs[k] > nt:
            return 0, 0, 1
            
    for k in range(1, nv + 1):
        for j in range(nr - 1):
            if Bs[j] == k:
                x[k-1] = B[j]
                break
                
    f = -B[nr-1] # Minimize
    return f, x, 0
import numpy as np

def dfpopt(fun, delfun, x0, alpha0, crit, kmax):
    """
    dfpopt: Quasi-Newton method (Davidon-Fletcher-Powell(DFP) algorithm)
    """
    x0 = np.array(x0, dtype=float)
    n = len(x0)
    nfmax = 10     # maximum function calls during a line search
    ni = 2         # indicate poor interval reductions before sectioning
    
    h = alpha0
    beta = 0.9
    xz = x0.copy()
    f = fun(xz)
    x1_ls = 0      # line search용 x1 (matlab 코드의 x1)
    f1_ls = f      # line search용 f1 (matlab 코드의 f1)
    
    k = 0
    nc = 0
    fold = f
    H = np.eye(n)
    
    while True:
        k += 1
        if k > kmax:
            print('Maximum possible iterations exceeded.')
            break
            
        x = xz.copy()
        delf = np.array(delfun(x), dtype=float)
        
        # 중단 조건: 그레이디언트의 노름이 기준치 이하일 때
        if np.linalg.norm(delf) <= crit:
            break
            
        # 탐색 방향 결정
        d = -H @ delf
        d_norm = np.linalg.norm(d)
        if d_norm > 1e-18:
            d = d / d_norm
        
        # Line Search (Quadratic Approximation)
        xs, fs = quadappx(fun, xz, x, d, x1_ls, f1_ls, h, nfmax, ni, crit)
        
        x1_ls = 0
        f1_ls = fs
        
        # 위치 업데이트
        xz_new = xz + xs * d
        delf0 = np.array(delfun(xz_new), dtype=float)
        
        # 그레이디언트 변화량 (gamma)
        gamma = delf0 - delf
        
        # DFP Matrix Update
        # 1. H * gamma * gamma' * H / (gamma' * H * gamma)
        d0_H = H @ gamma
        Qh = gamma @ d0_H
        if abs(Qh) > 1e-15:
            H = H - np.outer(d0_H, d0_H) / Qh
            
        # 2. delta * delta' / (delta' * gamma)
        delta = xs * d
        Pq = delta @ gamma
        if abs(Pq) > 1e-15:
            H = H + np.outer(delta, delta) / Pq
            
        xz = xz_new
        h = beta * xs
        fpr = fs
        
        # 중단 조건: 함수값 변화가 거의 없을 때
        if abs(fpr - fold) < crit:
            nc += 1
            if nc == ni:
                break
        else:
            nc = 0
            
        fold = fpr
        # 원본 코드의 nd (n+1) 리셋 로직은 nd 업데이트가 누락되어 있으나 구조 유지
        # if nd == n+1: H = np.eye(n)
        
    iter_count = k
    xopt = xz
    fopt = fs
    return xopt, fopt, iter_count

def quadappx(fun, xz, x, d, x1, f1, h, nfmax, ni, crit):
    """
    Quadratic approximation method for line search
    """
    tau = (np.sqrt(5) - 1) / 2
    x2 = x1 + h
    f2 = fun(xz + x2 * d)
    
    # 3점 패턴 찾기 (Bracketing)
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
            A = (x1 - x2) * (x1 - x3)
            B = (x2 - x1) * (x2 - x3)
            C = (x3 - x1) * (x3 - x2)
            # Quadratic fit point
            denom = (f1/A + f2/B + f3/C)
            if abs(denom) < 1e-20:
                x4 = (x1 + x3) / 2
            else:
                x4 = (f1*(x2+x3)/A + f2*(x1+x3)/B + f3*(x1+x2)/C) / denom / 2
        else:
            if x2 <= (x1 + x3) / 2:
                x4 = x2 + (1 - tau) * (x3 - x2)
            else:
                x4 = x3 - (1 - tau) * (x2 - x1)
            vs = 0
            
        # 최소 변화폭 설정
        dxs = sf * min(abs(x2 - x1), abs(x3 - x2))
        if abs(x4 - x1) < dxs:
            x4 = x1 + dxs
        elif abs(x4 - x3) < dxs:
            x4 = x3 - dxs
        elif abs(x4 - x2) < dxs:
            x4 = x2 - dxs if x2 > (x1 + x3) / 2 else x2 + dxs
            
        f4 = fun(xz + x4 * d)
        
        # 구간 업데이트
        if x4 > x2:
            if f4 >= f2:
                x3, f3 = x4, f4
            else:
                x1, f1, x2, f2 = x2, f2, x4, f4
        else:
            if f4 >= f2:
                x1, f1 = x4, f4
            else:
                x3, f3, x2, f2 = x2, f2, x4, f4
                
        snew = abs(x3 - x1)
        fmnew = (f1 + f2 + f3) / 3.0
        
        if abs(x3 - x1) <= crit:
            break
        if abs(fmnew - fmold) <= crit:
            wc += 1
            if wc == 2: break
        else:
            wc = 0
            
        if snew / sold > tau:
            vc += 1
            if vc == ni:
                vc = 0
                vs = 1
        else:
            vc = 0
            vs = 0
            
    return x2, f2

# --- 테스트 실행 (Rosenbrock function 예제) ---
if __name__ == "__main__":
    fun = lambda x: 100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
    delfun = lambda x: np.array([-400*x[0]*(x[1] - x[0]**2) + 2*(x[0] - 1), 
                                 200*(x[1] - x[0]**2)])
    
    x0 = [-1.2, 1.0]
    crit = 1e-6
    alpha0 = 1.0
    kmax = 1000
    
    xopt, fopt, iter_count = dfpopt(fun, delfun, x0, alpha0, crit, kmax)
    
    print(f"Optimal point: {xopt}")
    print(f"Function value at optimal point: {fopt}")
    print(f"Iterations: {iter_count}")
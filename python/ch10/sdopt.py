import numpy as np

def sdopt(fun, delfun, x0, alpha0, crit, kmax):
    """
    sdopt.m: Steepest descent method (최급강하법)
    
    Inputs:
    fun    : 목적 함수
    delfun : 목적 함수의 그레디언트(경사) 함수
    x0     : 시작점 (리스트 또는 배열)
    alpha0 : 초기 스텝 사이즈
    crit   : 중단 기준 (허용 오차)
    kmax   : 최대 반복 횟수
    
    Outputs:
    xopt   : 최적점
    fopt   : 최적점에서의 함수값
    iter   : 수행된 반복 횟수
    """
    xz = np.array(x0, dtype=float)
    nfmax = 10  # 라인 서치 중 최대 함수 호출 횟수
    ni = 2      # 구간 감소가 불충분할 때의 지표
    h = alpha0
    beta = 0.9
    
    f = fun(xz)
    x1 = 0.0
    f1 = f
    k = 0
    nc = 0
    fold = f
    fs = f
    
    while True:
        k += 1
        if k > kmax:
            print('Number of maximum iterations exceeded.')
            break
            
        x = xz.copy()
        delf = np.array(delfun(x))
        
        # 중단 조건: 그레디언트의 크기가 기준치보다 작을 때
        if np.linalg.norm(delf) <= crit:
            break
            
        # 최급강하 방향 벡터 (음의 그레디언트 방향)
        d = -delf
        norm_d = np.linalg.norm(d)
        if norm_d > 0:
            d = d / norm_d
            
        # 이차 근사법을 이용한 라인 서치 수행
        xs, fs = quadappx(fun, xz, x, d, x1, f1, h, nfmax, ni, crit)
        
        # 상태 업데이트
        x1 = 0.0
        f1 = fs
        xz = xz + xs * d
        h = beta * xs
        fpr = fs
        
        # 함수값 변화량에 따른 중단 조건
        if abs(fpr - fold) < crit:
            nc += 1
            if nc == ni:
                break
        else:
            nc = 0
            
        fold = fpr
        
    iter_count = k
    xopt = xz
    fopt = fs
    return xopt, fopt, iter_count

def quadappx(fun, xz, x, d, x1, f1, h, nfmax, ni, crit):
    """
    라인 서치를 위한 이차 근사법 (Quadratic approximation method)
    """
    tau = (np.sqrt(5) - 1) / 2
    x2 = x1 + h
    f2 = fun(xz + x2 * d)
    
    # 1. 초기 3점 패턴(3-point pattern) 찾기
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
                
    sf = 0.05  # 안전 계수 (0 < sf < 0.5)
    if (x1 >= x2 or x2 >= x3):
        # print('Incorrect interval.')
        return x2, f2
    if (f1 <= f2 or f2 >= f3):
        # print('Not 3-point pattern.')
        return x2, f2
        
    vs = 0
    vc = 0
    wc = 0
    j = 1
    
    # 2. 반복적인 이차 근사 최적화
    while j <= nfmax:
        sold = abs(x3 - x1)
        fmold = (f1 + f2 + f3) / 3.0
        
        if vs == 0:
            # 이차 보간법 식 적용
            A = (x1 - x2) * (x1 - x3)
            B = (x2 - x1) * (x2 - x3)
            C = (x3 - x1) * (x3 - x2)
            denom = (f1/A + f2/B + f3/C)
            if abs(denom) < 1e-20: # 분모 0 방지
                x4 = (x1 + x3) / 2
            else:
                x4 = (f1*(x2+x3)/A + f2*(x1+x3)/B + f3*(x1+x2)/C) / denom / 2
        else:
            # Golden section 형태의 안전 장치
            if x2 <= (x1 + x3) / 2:
                x4 = x2 + (1 - tau) * (x3 - x2)
            else:
                x4 = x3 - (1 - tau) * (x2 - x1)
            vs = 0
            
        # 인접 점들과 너무 겹치지 않게 보호 조치 (safeguard)
        dxs = sf * min(abs(x2 - x1), abs(x3 - x2))
        if abs(x4 - x1) < dxs:
            x4 = x1 + dxs
        elif abs(x4 - x3) < dxs:
            x4 = x3 - dxs
        elif abs(x4 - x2) < dxs:
            if x2 > (x1 + x3) / 2:
                x4 = x2 - dxs
            else:
                x4 = x2 + dxs
        
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
        j += 1
        
    return x2, f2

# --- 사용 예시 (Rosenbrock Function) ---
if __name__ == "__main__":
    # Example: Rosenbrock function
    def fun(x):
        return 100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
        
    def delfun(x):
        return np.array([
            -400*x[0]*(x[1] - x[0]**2) + 2*(x[0] - 1),
            200*(x[1] - x[0]**2)
        ])
    
    x0 = [-1.2, 1.0]
    crit = 1e-6
    alpha0 = 5.0
    kmax = 1000
    
    xopt, fopt, iters = sdopt(fun, delfun, x0, alpha0, crit, kmax)
    
    print(f"Optimal point: {xopt}")
    print(f"Function value: {fopt}")
    print(f"Iterations: {iters}")
import numpy as np

def saopt(fun, T, xl, xu, rp, rs, crit):
    """
    saopt.m: 시뮬레이티드 어닐링(SA) 방법을 이용한 최소화
    
    Inputs:
    fun  : 목적 함수
    T    : 초기 온도
    xl   : 하한값 배열 (lower bounds)
    xu   : 상한값 배열 (upper bounds)
    rp   : 온도 감소 계수 (temperature reduction factor)
    rs   : 스텝 감소 파라미터 (step reduction parameter)
    crit : 중단 기준 (stopping criterion)
    
    Outputs:
    xopt : 최적점
    fopt : 최적점에서의 함수값
    iter : 사이클 반복 횟수
    """
    # 초기 설정
    xl = np.array(xl, dtype=float)
    xu = np.array(xu, dtype=float)
    n = len(xu)
    
    sf = 2.0
    stp = 1.0
    np_iter = 10  # temperature step loop
    nc = 20       # search cycles
    nt = 2e4      # (원본 코드에 정의되어 있으나 루프에 사용되지 않음)
    ir = 16       # convergence check counter limit
    rcrit = 1e-10
    
    # 가능한 시작점 설정 (Random starting point)
    x = xl + np.random.rand(n) * (xu - xl)
    xs = x.copy()
    xmin = x.copy()
    
    f = fun(x)
    fmin = f
    fold = f
    citer = 0
    
    # 스텝 사이즈 및 수락률(acceptance ratio) 초기화
    stsize = np.full(n, stp)
    ar = np.ones(n)
    
    while True:
        for piter in range(np_iter): # 온도 단계 루프
            for ic in range(nc):    # 좌표 방향별 탐색 사이클
                for k in range(n):
                    # 새로운 후보 점 생성
                    xs[k] = x[k] + (2 * np.random.rand() - 1) * stsize[k]
                    
                    # 경계값 조건 확인
                    if (xs[k] < xl[k]) or (xs[k] > xu[k]):
                        xs[k] = xl[k] + np.random.rand() * (xu[k] - xl[k])
                    
                    fs = fun(xs)
                    
                    # 수락 여부 결정 (Metropolis criterion)
                    if fs <= f:
                        # 더 좋은 점이면 수락 및 최적값 업데이트
                        x[k] = xs[k]
                        f = fs
                        if fs < fmin:
                            xmin = xs.copy()
                            fmin = fs
                    else:
                        # 더 나쁜 점일 경우 확률적 수락
                        p = np.exp((f - fs) / T)
                        if np.random.rand() < p:
                            x[k] = xs[k]
                            f = fs
                        else:
                            # 거절 시 현재 좌표 복구 및 수락률 감소
                            xs[k] = x[k]
                            ar[k] = ar[k] - 1.0 / nc
            
            # 수락률이 약 50%가 되도록 스텝 사이즈 조정
            for j in range(n):
                if ar[j] > 0.6:
                    stsize[j] = stsize[j] * (1 + sf * (ar[j] - 0.6) / 0.4)
                elif ar[j] < 0.4:
                    stsize[j] = stsize[j] / (1 + sf * (0.4 - ar[j]) / 0.4)
                
                # 스텝 사이즈가 전체 범위를 넘지 않도록 제한
                if stsize[j] > (xu[j] - xl[j]):
                    stsize[j] = xu[j] - xl[j]
                
                ar[j] = 1.0 # 다음 사이클을 위해 수락률 초기화
        
        # 수렴 판정
        fcrit = crit + rcrit * abs(fmin)
        if (fmin <= fold) and (fold - fmin < fcrit):
            citer += 1
            if citer >= ir:
                break
        else:
            citer = 0
            
        # 온도 감소 및 상태 업데이트
        T = rp * T
        stp = rs * stp
        x = xmin.copy()
        f = fmin
        fold = f
        
    xopt = xmin
    fopt = fmin
    iter_out = citer
    
    return xopt, fopt, iter_out

# --- 사용 예시 (Example) ---
if __name__ == "__main__":
    import math
    
    # 예제 목적 함수 정의
    def f_example(x):
        # f = -cos(5*sqrt(sum(x-5).^2)) + 0.1*sum(x-5).^2
        dist_sq = np.sum((x - 5)**2)
        return -math.cos(5 * math.sqrt(dist_sq)) + 0.1 * dist_sq

    # 초기 파라미터 설정
    T_init = 100.0
    xl_val = np.array([-10.0, -10.0])
    xu_val = np.array([10.0, 10.0])
    rp_val = 0.8
    rs_val = 0.9
    crit_val = 1e-8

    # 함수 실행
    xopt, fopt, iter_count = saopt(f_example, T_init, xl_val, xu_val, rp_val, rs_val, crit_val)

    print(f"Optimal Point: {xopt}")
    print(f"Optimal Value: {fopt}")
    print(f"Iterations: {iter_count}")
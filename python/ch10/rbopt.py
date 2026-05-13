import numpy as np

def rbopt(fros, x0, crit):
    """
    rbopt.m: Rosenbrock’s method를 이용한 최소화 수치해석
    
    Inputs:
    fros  : 목적 함수 (콜백 함수)
    x0    : 시작점 (리스트 또는 넘파이 배열)
    crit  : 수렴 판정 기준
    
    Outputs:
    xopt  : 최적점
    fopt  : 최적점에서의 함수값
    istag : 수행된 스테이지(Stage) 수
    """
    x0 = np.array(x0, dtype=float)
    n = len(x0)  # 변수의 개수
    
    # 초기 설정
    mi = 7
    stsize = 1.0
    alpha = 3.0
    beta = 0.5
    xt = x0.copy()
    
    # 좌표축 방향의 초기 벡터 설정 (n x n 단위 행렬)
    v = np.eye(n)
    
    ft = fros(xt)
    istag = 0
    fold = ft
    
    while True:
        istag += 1
        
        # 각 스테이지 초기화
        stp = np.full(n, stsize) # 현재 단계의 스텝 크기
        stj = np.zeros(n)        # 각 방향으로의 성공적인 이동 거리 합계
        sfv = np.zeros((n, 2))   # 성공(column 0) 및 실패(column 1) 기록
        
        iter_idx = 0  # MATLAB의 iter (1 to n) -> 파이썬 0 to n-1
        icont = 0     # 성공과 실패가 모두 발생한 방향의 개수 추적
        idn = 0       # 연속 실패 횟수
        
        # 탐색 루프
        while True:
            # 탐색 방향 결정 (0부터 n-1까지 순환)
            x = xt + stp[iter_idx] * v[:, iter_idx]
            f = fros(x)
            
            if f < ft:
                # 성공 (Success)
                stj[iter_idx] += stp[iter_idx]
                xt = x.copy()
                ft = f
                idn = 0
                stp[iter_idx] *= alpha
                
                if sfv[iter_idx, 0] == 0:
                    sfv[iter_idx, 0] = 1
                    icont += 1
            else:
                # 실패 (Failure)
                stp[iter_idx] *= -beta
                
                if sfv[iter_idx, 1] == 0:
                    sfv[iter_idx, 1] = 1
                    icont += 1
                idn += 1
            
            # 모든 방향에서 최소 한 번씩 성공과 실패를 경험하면 스테이지 종료
            if icont == 2 * n:
                break
            
            # 비정상적인 무한 루프 방지
            if idn > 50:
                return xt, ft, istag
            
            # 다음 방향으로 인덱스 이동
            iter_idx = (iter_idx + 1) % n

        # 수렴 확인
        if istag > 1:
            if abs(ft - fold) < crit:
                print('Convergence achieved.')
                break
        
        fold = ft
        
        # Gram-Schmidt 과정을 위한 새로운 벡터 세트(u) 생성
        u = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                u[:, i] += stj[j] * v[:, j]
        
        # Gram-Schmidt 직교화 (Orthogonalization)
        for i in range(n):
            if i > 0:
                proj = np.zeros(n)
                for k in range(i):
                    c = np.dot(u[:, i], v[:, k])
                    proj += c * v[:, k]
                u[:, i] = u[:, i] - proj
            
            c = np.linalg.norm(u[:, i])
            
            # 첫 번째 벡터의 길이를 다음 스테이지의 기본 스텝 크기로 설정
            if i == 0:
                stsize = c
                if stsize < crit:
                    break
            
            # 직교성 상실 방지 및 정규화
            if c < 1e-15: # 매우 작은 값일 경우 좌표축 재설정
                v = np.eye(n)
                break
            else:
                v[:, i] = u[:, i] / c
                
    xopt = xt
    fopt = ft
    return xopt, fopt, istag

# 사용 예시 (Example)
if __name__ == "__main__":
    x0 = [-3.0, -1.0, 0.0, 1.0]
    crit = 1e-4
    
    # Powell's Quartic Function 예시
    def f(x):
        return (x[0] + 10*x[1])**2 + 5*(x[2] - x[3])**2 + (x[1] - 2*x[2])**4 + 10*(x[0] - x[3])**4

    xopt, fopt, istag = rbopt(f, x0, crit)
    
    print(f"Optimal Point: {xopt}")
    print(f"Function Value: {fopt}")
    print(f"Number of Stages: {istag}")
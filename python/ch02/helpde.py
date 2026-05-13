import numpy as np

def helpde(f, g, qx0, qxf, qy0, qyf, R, m, n, crit, kmax):
    """
    헬름홀츠 방정식 u_xx + u_yy + g(x,y)u = f(x,y)를 해결합니다.
    
    입력:
    f, g: (x, y)를 인자로 받는 함수
    qx0, qxf: x=x0, x=xf에서의 경계 조건 함수
    qy0, qyf: y=y0, y=yf에서의 경계 조건 함수
    R: [x0, xf, y0, yf] 영역 범위
    m, n: x축, y축 분할 수
    crit: 수렴 오차 허용치
    kmax: 최대 반복 횟수
    """
    x0, xf, y0, yf = R
    dx = (xf - x0) / m
    dy = (yf - y0) / n
    
    x = np.linspace(x0, xf, m + 1)
    y = np.linspace(y0, yf, n + 1)
    
    dx2, dy2 = dx**2, dy**2
    dxy2 = 2 * (dx2 + dy2)
    bx = dx2 / dxy2
    by = dy2 / dxy2
    bxy = (dx2 * dy2) / dxy2
    
    # u 행렬 초기화 (행: y좌표, 열: x좌표 구조가 시각화에 유리)
    # MATLAB의 u(i, j) 구조를 유지하기 위해 (m+1, n+1) 크기로 생성
    u = np.zeros((m + 1, n + 1))
    
    # 1. 경계 조건 설정
    # x = x0 (첫 번째 열), x = xf (마지막 열)
    for k in range(n + 1):
        u[0, k] = qx0(y[k])
        u[m, k] = qxf(y[k])
        
    # y = y0 (첫 번째 행), y = yf (마지막 행)
    for k in range(m + 1):
        u[k, 0] = qy0(x[k])
        u[k, n] = qyf(x[k])
        
    # 2. f(x,y)와 g(x,y) 값 미리 계산
    fv = np.array([[f(xi, yj) for yj in y] for xi in x])
    gv = np.array([[g(xi, yj) for yj in y] for xi in x])
    
    # 3. 반복법을 이용한 차분 방정식 풀이
    u0 = np.copy(u)
    for k in range(kmax):
        # 내부 격자점 업데이트 (경계 제외)
        for i in range(1, m):
            for j in range(1, n):
                u[i, j] = (by * (u[i+1, j] + u[i-1, j]) + 
                           bx * (u[i, j+1] + u[i, j-1]) + 
                           bxy * (gv[i, j] * u[i, j] - fv[i, j]))
        
        # 수렴 여부 확인
        if k > 0:
            diff = np.max(np.abs(u - u0))
            if diff < crit:
                print(f"Convergence reached at iteration {k}")
                break
        u0 = np.copy(u)
        
    return u, x, y
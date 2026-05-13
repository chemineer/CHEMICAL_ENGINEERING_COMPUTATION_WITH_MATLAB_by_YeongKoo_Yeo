import numpy as np

def parab1D(nx, nt, dx, dt, alpha, u0, bc, func=None, *args):
    """
    1차원 포물형 편미분 방정식(PDE)을 수치적으로 해석합니다.
    """
    # 초기화
    nx = int(nx)
    x = np.linspace(0, nx * dx, nx + 1)
    nt = int(nt)
    t = np.linspace(0, nt * dt, nt + 1)
    r = alpha * dt / dx**2
    
    u0 = np.array(u0).flatten()
    if len(u0) != nx + 1:
        raise ValueError('초기 조건 벡터의 길이가 맞지 않습니다.')
        
    bc = np.array(bc)
    if bc.shape[0] != 2:
        raise ValueError('경계 조건 행렬의 행 수는 2여야 합니다.')
    if bc.shape[1] < 2 or bc.shape[1] > 3:
        raise ValueError('경계 조건 행렬의 열 수는 2 또는 3이어야 합니다.')
    
    # bc가 2열인 경우 3열(0으로 채움) 추가
    if bc.shape[1] == 2:
        bc = np.hstack([bc, np.zeros((2, 1))])
        
    u = np.zeros((nx + 1, nt + 1))
    u[:, 0] = u0
    
    # 시간 반복
    for n in range(1, nt + 1):
        A = np.zeros((nx + 1, nx + 1))
        c = np.zeros(nx + 1)
        
        # 하단 x 경계 조건
        if bc[0, 0] == 1:
            A[0, 0] = 1
            c[0] = bc[0, 1]
        elif bc[0, 0] in [2, 3]:
            A[0, 0] = -3/(2*dx) - bc[0, 2]
            A[0, 1] = 2/dx
            A[0, 2] = -1/(2*dx)
            c[0] = bc[0, 1]
            
        # 내부 점
        for i in range(1, nx):
            A[i, i-1] = -r
            A[i, i] = 2*(1+r)
            A[i, i+1] = -r
            c[i] = r*u[i-1, n-1] + 2*(1-r)*u[i, n-1] + r*u[i+1, n-1]
            
            if func is not None:
                intercept = func(0, x[i], t[n], *args)
                slope = func(1, x[i], t[n], *args) - intercept
                A[i, i] -= dt * slope
                c[i] += dt * func(u[i, n-1], x[i], t[n-1], *args) + dt * intercept
                
        # 상단 x 경계 조건
        if bc[1, 0] == 1:
            A[nx, nx] = 1
            c[nx] = bc[1, 1]
        elif bc[1, 0] in [2, 3]:
            A[nx, nx] = 3/(2*dx) - bc[1, 2]
            A[nx, nx-1] = -2/dx
            A[nx, nx-2] = 1/(2*dx)
            c[nx] = bc[1, 1]
            
        u[:, n] = np.linalg.solve(A, c)
        
    return x, t, u
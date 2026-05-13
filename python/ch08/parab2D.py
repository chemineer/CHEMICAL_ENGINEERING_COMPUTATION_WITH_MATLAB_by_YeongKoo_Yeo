import numpy as np

def parab2D(nx, ny, nt, dx, dy, dt, alpha, u0, bc, func=None, *args):
    """
    2차원 포물형 편미분 방정식(PDE)을 명시적 방법(explicit method)으로 해석합니다.
    """
    nx = int(nx)
    ny = int(ny)
    x = np.linspace(0, nx * dx, nx + 1)
    y = np.linspace(0, ny * dy, ny + 1)
    
    # 안정성 조건 확인 및 dt 조정
    tmax = dt * nt if dt is not None else 0 
    if dt is None or dt > (dx**2 + dy**2) / (16 * alpha):
        dt = (dx**2 + dy**2) / (16 * alpha)
        nt = int(tmax / dt) + 1
        print(f"\ndt는 안정성을 위해 {dt:6.2e}로 조정되었습니다 (nt={nt:3d})")
        
    nt = int(nt)
    t = np.linspace(0, nt * dt, nt + 1)
    rx = alpha * dt / dx**2
    ry = alpha * dt / dy**2
    
    u0 = np.array(u0)
    if u0.shape != (nx + 1, ny + 1):
        raise ValueError('초기 조건 행렬의 크기가 맞지 않습니다.')
    
    bc = np.array(bc)
    if bc.shape[0] != 4:
        raise ValueError('경계 조건 행렬의 행 수는 4여야 합니다.')
        
    # bc가 2열인 경우 3열(0으로 채움) 추가
    if bc.shape[1] == 2:
        bc = np.hstack([bc, np.zeros((4, 1))])
        
    # 결과 배열 초기화
    u = np.zeros((nx + 1, ny + 1, nt + 1))
    u[:, :, 0] = u0
    
    # 시간 반복
    for n in range(nt):
        # 내부 점 계산
        for i in range(1, nx):
            for j in range(1, ny):
                u[i, j, n+1] = (rx * (u[i+1, j, n] + u[i-1, j, n]) + 
                                ry * (u[i, j+1, n] + u[i, j-1, n]) + 
                                (1 - 2*rx - 2*ry) * u[i, j, n])
                
                if func is not None:
                    u[i, j, n+1] += dt * func(u[i, j, n], x[i], y[j], t[n], *args)
        
        # 경계 조건 적용
        # Lower x
        if bc[0, 0] == 1:
            u[0, 1:ny, n+1] = bc[0, 1]
        else: # Neumann/Robbins
            u[0, 1:ny, n+1] = (-2*bc[0, 1]*dx + 4*u[1, 1:ny, n+1] - u[2, 1:ny, n+1]) / (2*bc[0, 2]*dx + 3)
            
        # Upper x
        if bc[1, 0] == 1:
            u[nx, 1:ny, n+1] = bc[1, 1]
        else:
            u[nx, 1:ny, n+1] = (-2*bc[1, 1]*dx - 4*u[nx-1, 1:ny, n+1] + u[nx-2, 1:ny, n+1]) / (2*bc[1, 2]*dx - 3)
            
        # Lower y
        if bc[2, 0] == 1:
            u[1:nx, 0, n+1] = bc[2, 1]
        else:
            u[1:nx, 0, n+1] = (-2*bc[2, 1]*dy + 4*u[1:nx, 1, n+1] - u[1:nx, 2, n+1]) / (2*bc[2, 2]*dy + 3)
            
        # Upper y
        if bc[3, 0] == 1:
            u[1:nx, ny, n+1] = bc[3, 1]
        else:
            u[1:nx, ny, n+1] = (-2*bc[3, 1]*dy - 4*u[1:nx, ny-1, n+1] + u[1:nx, ny-2, n+1]) / (2*bc[3, 2]*dy - 3)
            
        # 모서리 노드 평균화
        u[0, 0, n+1] = (u[0, 1, n+1] + u[1, 0, n+1]) / 2
        u[nx, 0, n+1] = (u[nx, 1, n+1] + u[nx-1, 0, n+1]) / 2
        u[0, ny, n+1] = (u[0, ny-1, n+1] + u[1, ny, n+1]) / 2
        u[nx, ny, n+1] = (u[nx, ny-1, n+1] + u[nx-1, ny, n+1]) / 2
        
    return x, y, t, u
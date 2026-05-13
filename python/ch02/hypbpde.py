import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def hypbpde(f1, f2, g0, g1, xspan, tspan, nx, nt, alpa):
    """
    1차원 쌍곡선형 PDE (파동 방정식: u_tt = alpa * u_xx) 풀이
    """
    x0, xf = xspan
    t0, tf = tspan
    
    dx = (xf - x0) / nx
    dt = (tf - t0) / nt
    
    # x는 (nx+1, 1) 형태, t는 (nt+1,) 형태
    x = np.linspace(x0, xf, nx + 1)
    t = np.linspace(t0, tf, nt + 1)
    
    # 파라미터 q 계산
    q = alpa * (dt / dx)**2
    q1 = q / 2
    q2 = 2 * (1 - q)
    
    # 결과 행렬 u 초기화 (행: x, 열: t)
    u = np.zeros((nx + 1, nt + 1))
    
    # 1. 초기 조건 (t = 0)
    u[:, 0] = [f1(xi) for xi in x]
    
    # 2. 경계 조건 (x = x0, x = xf)
    for k in range(nt + 1):
        u[0, k] = g0(t[k])
        u[nx, k] = g1(t[k])
        
    # 3. 첫 번째 시간 단계 계산 (Central difference for u_t 초기조건 반영)
    # MATLAB: u(2:nx, 2)
    u[1:nx, 1] = (q1 * u[0:nx-1, 0] + 
                  (1 - q) * u[1:nx, 0] + 
                  q1 * u[2:nx+1, 0] + 
                  dt * np.array([f2(xi) for xi in x[1:nx]]))
    
    # 4. 나머지 시간 단계 루프 (k=3부터 시작)
    # Python range(2, nt+1)은 k=2부터 nt까지 반복 (MATLAB의 3:nt+1에 대응)
    for k in range(2, nt + 1):
        u[1:nx, k] = (q * u[0:nx-1, k-1] + 
                      q2 * u[1:nx, k-1] + 
                      q * u[2:nx+1, k-1] - 
                      u[1:nx, k-2])
        
    # 5. 결과 시각화 (surf)
    T, X = np.meshgrid(t, x)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(T, X, u, cmap='gray', edgecolor='none')
    
    ax.set_xlabel('t')
    ax.set_ylabel('x')
    ax.set_zlabel('u(t,x)')
    plt.colorbar(surf)
    plt.show()
    
    return u, q

# --- 사용 예시 ---
# f1: 초기 변위, f2: 초기 속도, g0/g1: 경계 고정
# result_u, q_val = hypbpde(lambda x: np.sin(np.pi*x), lambda x: 0, 
#                           lambda t: 0, lambda t: 0, 
#                           [0, 1], [0, 2], 20, 50, 1)
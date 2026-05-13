import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def parapde(f, g0, g1, tf, nx, nt, alpa):
    """
    명시적 방법을 이용한 1차원 포물선형 PDE (u_t = alpa * u_xx) 풀이
    """
    # 1. 그리드 설정
    h = 1 / nx
    d = tf / nt
    r = alpa * d / h**2
    
    # 안정성 조건 체크
    if r > 0.5:
        print(f"Warning: r = {r:.4f}가 0.5보다 큽니다. 수치적 불안정성이 발생할 수 있습니다.")
    
    x = np.linspace(0, 1, nx + 1)
    t = np.linspace(0, tf, nt + 1)
    
    # 결과 행렬 u 초기화 (행: x, 열: t)
    u = np.zeros((nx + 1, nt + 1))
    
    # 2. 초기 조건 설정: u(x, 0) = f(x)
    u[:, 0] = f(x)
    
    # 3. 경계 조건 설정: u(0, t) = g0(t), u(1, t) = g1(t)
    u[0, :] = g0(t)
    u[nx, :] = g1(t)
    
    # 4. 시간 전진 루프 (Explicit Method)
    # MATLAB: u(2:nx, k+1) = r*u(1:nx-1, k) + (1-2r)*u(2:nx, k) + r*u(3:nx+1, k)
    for k in range(nt):
        u[1:nx, k+1] = (r * u[0:nx-1, k] + 
                        (1 - 2*r) * u[1:nx, k] + 
                        r * u[2:nx+1, k])
    
    # 5. 결과 시각화 (MATLAB의 surf와 동일하게 u를 전치하여 시각화)
    U_display = u.T  # (nt+1, nx+1) 형태로 변환
    X, T = np.meshgrid(x, t)
    
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, T, U_display, cmap='gray', edgecolor='none')
    
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u(x,t)')
    plt.colorbar(surf)
    plt.show()
    
    return u

# --- 사용 예시 ---
# f = lambda x: np.sin(np.pi * x)  # 초기 온도 분포
# g0 = lambda t: np.zeros_like(t)  # 왼쪽 끝 고정 온도 0
# g1 = lambda t: np.zeros_like(t)  # 오른쪽 끝 고정 온도 0
# sol = parapde(f, g0, g1, 0.1, 20, 100, 1)
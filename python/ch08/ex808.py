import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 파라미터 설정
R = 0.03
tmax = 0.5
m = 1
N_x = 40  # r 방향 분할 수
N_t = 20  # L (t) 방향 분할 수

x = np.linspace(0, R, N_x)
dx = x[1] - x[0]
t_eval = np.linspace(0, tmax, N_t)

# 물성치 및 조건
Cp = 4.2
rho = 1e3
vk = 5.0
k = 0.1
Q = 100.0
T0 = 500.0

# 2. 미분방정식 정의 (Method of Lines)
def heat_equation(t, U):
    # U는 r=0 부터 r=R 직전까지의 내부 노드 온도 (크기: N_x - 1)
    dUdt = np.zeros_like(U)
    
    # [Zero Division 해결 2] r = R 경계조건 처리
    # c(R) = 0 이므로 ODE 통합에서 제외하고, 경계 조건(qr=-k, pr=-Q)을 이용해 대수적으로 계산합니다.
    # 2차 후진 차분(2nd-order backward difference)을 사용하여 u_R 계산: du/dx = -Q/k
    u_R = (4 * U[-1] - U[-2] - 2 * dx * Q / k) / 3
    
    # 전체 배열 결합 (인덱싱을 편하게 하기 위함)
    u = np.append(U, u_R)
    
    # [Zero Division 해결 1] r = 0 (i = 0) 처리
    # 원통형 좌표계의 1/r 항으로 인한 0 나누기를 방지하기 위해 로피탈의 정리(L'Hopital's rule) 적용
    c_0 = rho * Cp * vk * (1 - (x[0] / R)**2)
    # 1/r * d/dr(r * du/dr) -> 2 * d^2u/dr^2 로 수렴. 
    # 대칭성(u[-1] = u[1])을 이용하여 차분.
    dUdt[0] = (1 / c_0) * 4 * (u[1] - u[0]) / (dx**2)
    
    # 0 < r < R 내부 노드 (i = 1 ~ N_x - 2) 처리
    for i in range(1, len(U)):
        xi = x[i]
        c_i = rho * Cp * vk * (1 - (xi / R)**2)
        
        # 중심 차분법(Central difference)을 이용한 플럭스 계산
        x_plus = xi + dx / 2
        x_minus = xi - dx / 2
        
        flux_out = x_plus * (u[i+1] - u[i]) / dx
        flux_in = x_minus * (u[i] - u[i-1]) / dx
        
        dUdt[i] = (1 / (c_i * xi)) * (flux_out - flux_in) / dx
        
    return dUdt

# 3. 초기 조건 설정 (내부 노드만)
U0 = np.full(N_x - 1, T0)

# 4. 상미분방정식(ODE) 풀이
# 확산 방정식은 Stiff한 특성이 있으므로 'BDF' (Backward Differentiation Formula) 메서드 사용
sol = solve_ivp(heat_equation, [0, tmax], U0, t_eval=t_eval, method='BDF')

# 5. 풀이 결과에 r=R 경계값 다시 추가하여 전체 해 구성
u_full = np.zeros((N_t, N_x))
for j in range(N_t):
    U_t = sol.y[:, j]
    u_R = (4 * U_t[-1] - U_t[-2] - 2 * dx * Q / k) / 3
    u_full[j, :] = np.append(U_t, u_R)

# 6. 시각화 (MATLAB의 subplot, mesh, plot 대응)
X, T_grid = np.meshgrid(x, t_eval)

fig = plt.figure(figsize=(12, 5))

# 3D Surface Plot
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax1.plot_surface(X, T_grid, u_full, cmap='jet', rstride=1, cstride=1, antialiased=True)
ax1.set_xlabel('r (m)')
ax1.set_ylabel('L (m)')
ax1.set_zlabel('T (K)')
ax1.set_title('Surface Plot')

# 2D Line Plot
ax2 = fig.add_subplot(1, 2, 2)
for j in range(N_t):
    ax2.plot(x, u_full[j, :], 'k')
ax2.set_xlabel('r (m)')
ax2.set_ylabel('T (K)')
ax2.grid(True)
ax2.set_title('T vs r(m)')

plt.tight_layout()
plt.show()
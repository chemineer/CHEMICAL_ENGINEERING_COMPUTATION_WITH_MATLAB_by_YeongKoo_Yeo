import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D

# 1. 파라미터 및 격자 설정
alpha = 4.8e-7
L = 1.0           # x의 범위 (0에서 1)
nx = 20           # 공간 격자 수
x = np.linspace(0, L, nx)
dx = x[1] - x[2]

t_span = (0, 21600)  # 시작 시간, 종료 시간
t_eval = np.linspace(0, 21600, 54) # 결과 출력 시간대

# 2. 초기 조건 설정 (Initial Condition)
# 모든 지점의 초기 온도 T = 90
T0 = np.full(nx, 90.0)

# 3. PDE 정의 (공간에 대한 이산화)
def heat_equation(t, T):
    dTdt = np.zeros_like(T)
    
    # 경계 조건 (Boundary Conditions): pl=ul-15, pr=ur-15 (Dirichlet 조건)
    # 양쪽 끝 온도를 15도로 고정
    T[0] = 15
    T[-1] = 15
    
    # 내부 점들에 대한 중앙 차분 (Central Difference)
    # dT/dt = alpha * d^2T/dx^2
    for i in range(1, nx - 1):
        d2Tdx2 = (T[i+1] - 2*T[i] + T[i-1]) / (dx**2)
        dTdt[i] = alpha * d2Tdx2
        
    # 경계점의 변화율은 0으로 유지 (온도가 고정되어 있으므로)
    dTdt[0] = 0
    dTdt[-1] = 0
    
    return dTdt

# 4. 수치 해석 실행 (Solver)
sol = solve_ivp(heat_equation, t_span, T0, t_eval=t_eval, method='RK45')

# 5. 결과 시각화 (3D Surface Plot)
X, T_grid = np.meshgrid(x, sol.t)
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# sol.y.T는 시간(row) x 공간(col) 배열입니다.
surf = ax.plot_surface(X, T_grid, sol.y.T, cmap='viridis')

ax.set_xlabel('x')
ax.set_ylabel('t(sec)')
ax.set_zlabel('T(deg C)')
ax.set_title('Temperature Profile')
plt.show()
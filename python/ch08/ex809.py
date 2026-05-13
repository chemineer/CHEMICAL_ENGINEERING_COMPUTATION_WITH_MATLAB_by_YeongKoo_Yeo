import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pflowht import pflowht  # 작성하신 pflowht.py를 임포트

# 1. 데이터 및 파라미터 설정 (Data and parameters)
L = 2
R = 0.05
n = 20
h = R / n
r = np.arange(0, n + 1) * h
alpa = 1e-4
Ti = 300
Tb = 400
vmax = 0.5

# 속도 분포 계산 (v(R)=0 이 될 경우 Zero Division 발생 가능하므로 주의)
v = vmax * (1 - (r / R)**2)

# pflowht.py의 pf 인자를 충족하기 위한 객체 생성
class Parameter:
    def __init__(self, r, v, n, h, alpa, Tb):
        self.r = r
        self.v = v
        self.n = n
        self.h = h
        self.alpa = alpa
        self.Tb = Tb

pf = Parameter(r, v, n, h, alpa, Tb)

# 2. PDE 해결 (Solve PDE)
# 초기 조건: 모든 노드의 온도를 Ti로 설정 (size: n)
T0 = np.ones(n) * Ti

# MATLAB의 ode15s에 대응하는 'BDF' 또는 'Radau' 메서드 사용 (Stiff 문제)
# rtol과 atol은 MATLAB의 odeset 설정값을 반영
sol = solve_ivp(
    lambda x, T: pflowht(x, T, pf), 
    [0, L], 
    T0, 
    method='BDF', 
    rtol=1e-12, 
    atol=1e-10
)

# 3. 결과 정리 (Display results)
x_sol = sol.t  # x 방향 (Time step in ODE)
T_sol = sol.y.T  # 각 위치 x에서의 온도 분포 (Internal nodes)

# 벽면 온도(Tb) 열 추가 (MATLAB의 T = [T, Tw] 과정)
na = len(x_sol)
Tw = np.ones((na, 1)) * Tb
T_final = np.hstack((T_sol, Tw))

# 4. 시각화 (Subplots)
fig = plt.figure(figsize=(12, 5))

# Subplot 1: Mesh Plot (3D Surface)
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
R_grid, X_grid = np.meshgrid(r, x_sol)
surf = ax1.plot_surface(R_grid, X_grid, T_final, cmap='gray', edgecolor='none')
ax1.set_xlabel('r(m)')
ax1.set_ylabel('x(m)')
ax1.set_zlabel('T(K)')
ax1.set_title('Temperature Profile (3D)')

# Subplot 2: Contour Plot
ax2 = fig.add_subplot(1, 2, 2)
# MATLAB의 contourf(x, r, T')와 동일하게 축 맞춤
cp = ax2.contourf(x_sol, r, T_final.T, cmap='gray')
fig.colorbar(cp, ax=ax2)
ax2.set_xlabel('x(m)')
ax2.set_ylabel('r(m)')
ax2.set_title('Temperature Contour')

plt.tight_layout()
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 공간 및 시간 설정
nx = 100
x = np.linspace(0, 1, nx)
t = np.linspace(0, 1, 60)
dx = x[1] - x[0]

# 2. PDE 정의: du/dt = d^2u/dx^2 + s
# MATLAB의 fleqn: c=1, f=DuDx, s=4 => du/dt = d^2u/dx^2 + 4
def pde_system(t, u):
    dudt = np.zeros_like(u)
    # 중앙 차분법 (Central Difference)을 통한 2계 미분
    d2udx2 = np.gradient(np.gradient(u, dx), dx)
    dudt = d2udx2 + 4
    
    # 경계 조건 적용
    dudt[0] = 0  # pl=0, ql=1 -> du/dx = 0
    dudt[-1] = 0 # pr=ur, qr=0 -> u(1) = 0
    return dudt

# 3. 초기 조건 및 경계 조건
u0 = np.zeros(nx)

# 4. PDE 풀이
sol = solve_ivp(pde_system, [0, 1], u0, t_eval=t)
u = sol.y.T # (time, space)

# 5. 시각화
fig = plt.figure(figsize=(12, 5))

# 좌측: Surf plot
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
X, T = np.meshgrid(x, t)
ax1.plot_surface(X, T, u, cmap='gray')
ax1.set_xlabel(r'$\xi$'); ax1.set_ylabel(r'$\tau$'); ax1.set_zlabel(r'$\phi$')

# 우측: Plot tau vs xi
ax2 = fig.add_subplot(1, 2, 2)
for k in range(len(t)):
    ax2.plot(x, u[k, :], 'k')
ax2.set_xlabel(r'$\xi$'); ax2.set_ylabel(r'$\tau$')

plt.tight_layout()
plt.show()
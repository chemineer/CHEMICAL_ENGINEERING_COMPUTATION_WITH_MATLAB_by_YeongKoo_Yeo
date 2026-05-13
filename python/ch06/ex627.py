import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pfrdiff import pfrdiff  # pfrdiff.py 임포트

# 1. 데이터 및 파라미터 설정
n = 50
pf = {'n': n, 'Pe': 1, 'Da': 2}
h = 1.0 / n

# 2. 초기화
Z0 = np.ones(n)         # 초기 농도 분포 (phi = 1)
tspan = [0, 1]          # 시간 구간 (tau)

# 3. ODE 세트 풀이 (ode45에 해당)
sol = solve_ivp(
    pfrdiff, 
    tspan, 
    Z0, 
    args=(pf,), 
    method='RK45', 
    t_eval=np.linspace(tspan[0], tspan[1], 100) # 그래프를 위한 시간 샘플링
)

t = sol.t               # 시간 벡터
C = sol.y.T             # 농도 행렬 (시간 x 위치)

# 4. 결과 정리[cite: 4]
# 정상 상태 농도 (마지막 시간 타임스텝)
Cs = np.insert(C[-1, :], 0, 1.0) 
x = np.linspace(h, 1.0, n)

# 5. 시각화 (MATLAB의 mesh에 해당)[cite: 4]
X, T = np.meshgrid(x, t)
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# 3D 표면 그래프 작성
surf = ax.plot_surface(X, T, C, cmap='viridis', edgecolor='none')

ax.set_xlabel(r'$\xi$ (Position)')
ax.set_ylabel(r'$\tau$ (Time)')
ax.set_zlabel(r'$\phi (\tau, \xi)$ (Concentration)')
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.title('PFR with Axial Diffusion (MoL)')
plt.show()

# 최종 농도 출력
print(f"Steady-state concentration at exit: {C[-1, -1]:.4f}")
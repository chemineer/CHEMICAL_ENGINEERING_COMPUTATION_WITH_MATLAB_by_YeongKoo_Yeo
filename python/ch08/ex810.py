import numpy as np
import matplotlib.pyplot as plt
from parabDbc import parabDbc  # 제공된 parabDbc.py 임포트

# 1. 데이터 및 파라미터 설정 (MATLAB ex810.m 로직 반영)
alpha = 4.52e-7        # 열확산 계수
nx = 25                # 공간 분할 수
dx = 0.5 / nx          # 공간 간격
dt = 300               # 시간 간격 (초)
nt = int(18000 / dt)   # 시간 단계 수

# 초기 조건 및 경계 조건
u0 = 100 * np.ones(nx + 1)  # 초기 온도 100
bci = 18                    # x=0에서의 경계 조건
bcf = 18                    # x=L에서의 경계 조건

# 2. PDE 해석 수행
# u: 온도 결과 행렬, r: 확산 계수 (alpha * dt / dx^2)
u, r_val = parabDbc(nx, nt, dx, dt, alpha, u0, bci, bcf)

print(f"Stability parameter r: {r_val:.4f}")

# 3. 결과 시각화 (MATLAB의 surf 및 view 설정 재현)
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# 그리드 생성
x_coords = np.arange(nx + 1)
t_coords = np.arange(nt + 1)
X, T = np.meshgrid(x_coords, t_coords)

# Surface plot 생성
# MATLAB의 colormap과 유사한 효과를 위해 'viridis' 혹은 'jet' 사용
surf = ax.plot_surface(X, T, u, cmap='jet', edgecolor='none', alpha=0.9)

# 축 라벨 및 범위 설정
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('T')
ax.set_zlim(0, 110)
ax.set_title('1D Heat Equation (Dirichlet Boundary Conditions)')

# MATLAB의 view([-217 30]) 설정 반영
# Matplotlib의 azimuth는 MATLAB과 기준이 다르므로 조정이 필요함
ax.view_init(elev=30, azim=-217)

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.tight_layout()
plt.show()
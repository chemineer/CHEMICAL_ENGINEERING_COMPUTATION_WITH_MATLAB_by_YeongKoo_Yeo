import numpy as np
import matplotlib.pyplot as plt
from ellipticPDE import ellipticPDE  # 제공된 ellipticPDE.py 임포트

# 1. 데이터 및 파라미터 설정 (MATLAB ex811.m 로직 반영)
distx = 1.0  # x 방향 플레이트 길이 (m)
disty = 1.0  # y 방향 플레이트 너비 (m)
ndx = 20     # x 방향 분할 수
ndy = 20     # y 방향 분할 수
rhf = -100e3 / 16  # 방정식의 우변 (Poisson 상수)

# 2. 경계 조건 설정 (Boundary conditions)
# 각 행: [유형, 상수(beta), 계수(gamma)]
# 유형 3: Robbins (혼합) 경계 조건
bc = np.array([
    [3, -5 * 25,  5], # Lower x boundary
    [3,  5 * 25, -5], # Upper x boundary
    [3, -5 * 25,  5], # Lower y boundary
    [3,  5 * 25, -5]  # Upper y boundary
])

# 3. PDE 해석 수행
# x, y: 좌표 그리드 벡터, T: 계산된 온도 행렬
dx = distx / ndx
dy = disty / ndy
x, y, T = ellipticPDE(ndx, ndy, dx, dy, bc, rhf)

# 4. 결과 시각화 (MATLAB의 surf 및 view 설정 재현)
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# 시각화를 위한 그리드 생성 (y, x 순서 주의)
X_grid, Y_grid = np.meshgrid(x, y)

# Surface plot 생성 (MATLAB의 colormap('jet') 반영)
surf = ax.plot_surface(X_grid, Y_grid, T, cmap='jet', edgecolor='none')

# 축 라벨 및 설정
ax.set_xlabel('x(m)')
ax.set_ylabel('y(m)')
ax.set_zlabel('T(deg C)')
ax.set_title('2D Elliptic PDE - Temperature Distribution')

# MATLAB의 view(135, 45) 설정 반영
ax.view_init(elev=45, azim=135)

# 컬러바 추가
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.tight_layout()
plt.show()
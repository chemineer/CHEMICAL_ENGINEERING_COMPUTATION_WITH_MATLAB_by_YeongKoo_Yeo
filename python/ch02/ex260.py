import numpy as np
import matplotlib.pyplot as plt
from helpde import helpde

# 1. 영역 및 기본 함수 설정
R = [0, 4, 0, 4]  # [x0, xf, y0, yf]
f = lambda x, y: 0
g = lambda x, y: 0

# 2. 경계 조건 함수 정의 (MATLAB의 익명 함수 대응)
qx0 = lambda y: np.exp(y) - np.cos(y)
qxf = lambda y: np.exp(y) * np.cos(4) - np.exp(4) * np.cos(y)
qy0 = lambda x: np.cos(x) - np.exp(x)
qyf = lambda x: np.exp(4) * np.cos(x) - np.exp(x) * np.cos(4)

# 3. 파라미터 설정
m = 40
n = 40
crit = 1e-6
kmax = 1000

# 4. helpde 함수 호출 (수치 해 계산)
u, x, y = helpde(f, g, qx0, qxf, qy0, qyf, R, m, n, crit, kmax)

# 5. 결과 시각화 (MATLAB의 mesh 함수 대응)
X, Y = np.meshgrid(x, y)
# helpde.py의 u는 u[x_idx, y_idx] 구조이므로 시각화를 위해 전치(Transpose)가 필요할 수 있음
# meshgrid와 일치시키기 위해 u.T 사용
Z = u.T 

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# mesh/surface 그래프 생성
surf = ax.plot_wireframe(X, Y, Z, color='gray', linewidth=0.5)

# 축 라벨 및 범위 설정
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u(x,y)')
ax.set_zlim(-100, 100)
ax.set_title('Steady-state Temperature Distribution')

plt.show()
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# 1. 미분 방정식 정의 (zode)
# dz/dx = [z1, -abs(z0)]
def zode(x, z):
    return np.vstack((z[1], -np.abs(z[0])))

# 2. 경계 조건 정의 (zobc)
# res = [za(0), zb(0) + 2]
def zobc(za, zb):
    return np.array([za[0], zb[0] + 2])

# 3. 초기 그리드 및 추정값 설정 (bvpinit)
# x는 0부터 4까지 10개의 구간, 초기 추정값은 [1, 0]
x_init = np.linspace(0, 4, 10)
z_init = np.zeros((2, x_init.size))
z_init[0, :] = 1
z_init[1, :] = 0

# 4. BVP 문제 해결 (bvp4c)
sol = solve_bvp(zode, zobc, x_init, z_init)

# 5. 결과 시각화 (deval & plot)
if sol.success:
    x_plot = np.linspace(0, 4, 100)
    y_plot = sol.sol(x_plot)[0]  # z(1,:) 에 해당

    plt.figure(figsize=(8, 5))
    plt.plot(x_plot, y_plot)
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('BVP Solution using solve_bvp')
    plt.show()
else:
    print("해를 찾지 못했습니다:", sol.message)
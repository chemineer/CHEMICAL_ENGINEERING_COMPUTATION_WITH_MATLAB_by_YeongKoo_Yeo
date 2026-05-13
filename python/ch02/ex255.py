import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# 1. 미분 방정식 시스템 정의 (y'' = 6y/x^2)
# z[0] = y, z[1] = y'
def bcprob(x, z):
    # dz/dx = [z[1], 6*z[0]/x^2]
    return np.array([z[1], 6 * z[0] / x**2])

# 2. 경계 조건 정의 (Boundary Conditions)
# za: 왼쪽 경계(x=1)에서의 값, zb: 오른쪽 경계(x=2)에서의 값
def bcval(za, zb):
    # y(1) - 1 = 0, y(2) - 1 = 0 (즉, y(1)=1, y(2)=1)
    return np.array([za[0] - 1, zb[0] - 1])

# 3. 초기 추측값(Initial Guess) 설정
# x: 1부터 2 사이를 10개 구간으로 나눔
# y: 모든 지점에서 y=1, y'=0으로 가정
x_init = np.linspace(1, 2, 10)
y_init = np.zeros((2, x_init.size))
y_init[0, :] = 1  # y의 초기 추측값
y_init[1, :] = 0  # y'의 초기 추측값

# 4. BVP 풀이
res = solve_bvp(bcprob, bcval, x_init, y_init)

# 5. 결과 시각화
x_plot = np.linspace(1, 2, 100)
y_plot = res.sol(x_plot)[0] # 첫 번째 행이 y(x) 값

plt.figure(figsize=(7, 5))
plt.plot(x_plot, y_plot)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution of Boundary Value Problem (BVP)')
plt.show()
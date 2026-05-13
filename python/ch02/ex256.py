import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp

# 1. 미분 방정식 시스템 정의
def metfun(x, y):
    """
    y[0]: 온도 T
    y[1]: 온도의 기울기 dT/dx
    """
    h = 50      # 열전달 계수
    D = 0.04    # 막대 직경
    k = 390     # 열전도도
    Ta = 25     # 주변 온도
    
    dydx = [y[1], 4*h*(y[0] - Ta) / (D * k)]
    return np.array(dydx)

# 2. 경계 조건 정의
def metbc(ya, yb):
    """
    ya: x=0 (왼쪽 끝)에서의 상태
    yb: x=1 (오른쪽 끝)에서의 상태
    """
    # ya[0] - 100 = 0  => T(0) = 100
    # yb[0] - 0 = 0    => T(1) = 0
    return np.array([ya[0] - 100, yb[0]])

# 3. 초기 추측값(Initial Guess) 및 구간 설정
# [0, 1] 구간을 20개의 하위 구간으로 나눔
x_init = np.linspace(0, 1, 20)
y_init = np.zeros((2, x_init.size))
y_init[0, :] = 100  # 온도 T에 대한 초기 추측값
y_init[1, :] = 0    # 기울기 dT/dx에 대한 초기 추측값

# 4. BVP 풀이
sol = solve_bvp(metfun, metbc, x_init, y_init)

# 5. 결과 시각화
x_plot = np.linspace(0, 1, 100)
y_plot = sol.sol(x_plot)[0]  # T(x) 값 추출

plt.figure(figsize=(7, 5))
plt.plot(x_plot, y_plot, label='Temperature Profile')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('T(degC)')
plt.title('Metal Rod Temperature Distribution (BVP)')
plt.legend()
plt.show()
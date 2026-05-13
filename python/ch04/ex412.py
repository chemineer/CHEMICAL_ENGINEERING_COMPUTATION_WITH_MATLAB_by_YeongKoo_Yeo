import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from bzacfun import bzacfun

# 1. 데이터 설정
x1 = np.array([0.0, 0.0069, 0.1565, 0.3396, 0.4666, 0.6004, 0.7021, 0.8286, 0.8862, 0.9165, 0.9561, 0.9840, 1.0])
P = np.array([57.52, 58.2, 126.00, 175.30, 189.50, 224.30, 236.00, 250.20, 259.00, 261.11, 264.45, 266.53, 271.00])

# 2. 다항식 회귀 및 미분
c = np.polyfit(x1, P, 4)
dc = np.polyder(c)

# 3. ODE 풀이
x10 = 1e-5
# [가장 중요한 수정점] MATLAB의 dc(end)/c(end)와 동일하게 파이썬 배열의 마지막 요소(-1)를 호출
y10 = x10 * dc[-1] / c[-1]

# MATLAB ode45와 동일하게 RK45 사용, 자동 스텝 사이즈를 위해 t_eval 생략
sol = solve_ivp(bzacfun, [x10, 1.0], [y10], method='RK45')

x = sol.t
y = sol.y[0]

# 4. 활동도 계수 계산
P1v = P[-1]
P2v = P[0]
Px = np.polyval(c, x)

# 파이썬에서는 분모가 0이 될 때(x=1) 경고가 뜨지만, MATLAB처럼 NaN으로 처리하고 넘어갑니다.
gam1 = (y * Px) / (x * P1v)
gam2 = ((1 - y) * Px) / ((1 - x) * P2v)

# 5. 결과 그래프 출력
plt.figure(figsize=(10, 5))

# x1-y1 그래프
plt.subplot(1, 2, 1)
plt.plot(x, y)
plt.xlabel('x_1')
plt.ylabel('y_1')

# 활동도 계수 그래프
plt.subplot(1, 2, 2)
plt.plot(x, gam1, label='$\gamma_1$')
plt.plot(x, gam2, '.-', label='$\gamma_2$')
plt.xlabel('x_1')
plt.ylabel('$\gamma$')
plt.legend(loc='best')
plt.axis([0, 1, 0, 4])

plt.tight_layout()
plt.show()
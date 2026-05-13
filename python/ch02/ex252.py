import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 파라미터 및 초기값 설정
a = [13.2, 0.95, 1.76]  # 계수 [a1, a2, a3]
x0 = [0.028, 0]         # 초기 조건 [x1(0), x2(0)]
t_span = (0, 1)         # 시간 범위
t_eval = np.linspace(0, 1, 100)  # 그래프를 매끄럽게 그리기 위한 시간 지점들

# 2. 미분 방정식 정의 (dx/dt)
def penrxn(t, x):
    x1, x2 = x
    # dx1/dt = a1*x1 - (a1/a2)*x1^2 (로지스틱 성장 모델)
    # dx2/dt = a3*x1 (생성물 형성 모델)
    dx1dt = a[0] * x1 - (a[0] / a[1]) * x1**2
    dx2dt = a[2] * x1
    return [dx1dt, dx2dt]

# 3. ODE 풀이 (MATLAB의 ode45와 유사)
sol = solve_ivp(penrxn, t_span, x0, t_eval=t_eval)

# 4. 결과 시각화
plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0], label='x1 (Biomass)', linewidth=1.5)
plt.plot(sol.t, sol.y[1], '--', label='x2 (Penicillin)', linewidth=1.5)

plt.grid(True)
plt.xlabel('t')
plt.ylabel('x1 and x2')
plt.title('Penicillin Production Reaction Profile')
plt.legend(loc='best')
plt.show()
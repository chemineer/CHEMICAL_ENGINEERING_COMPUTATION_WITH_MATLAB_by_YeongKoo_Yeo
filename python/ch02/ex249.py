import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 1. van der Pol 방정식 정의
# y[0] = y1, y[1] = y2 (즉, y1의 미분값)
def vdp_system(t, y, mu):
    dy1_dt = y[1]
    dy2_dt = -y[0] + mu * (1 - y[0]**2) * y[1]
    return [dy1_dt, dy2_dt]

# 2. 파라미터 및 초기값 설정
mu = 1
tspan = (0, 25)
y0 = [1, 1]

# 3. ODE 해결 (MATLAB의 ode45와 동일한 RK45 방식)
# args=(mu,)를 통해 함수에 파라미터를 전달합니다.
sol = solve_ivp(vdp_system, tspan, y0, args=(mu,), t_eval=np.linspace(0, 25, 500))

# 4. 결과 시각화
plt.figure(figsize=(10, 5))
plt.plot(sol.t, sol.y[0], label='y_1')           # 실선
plt.plot(sol.t, sol.y[1], ':', label='y_2')      # 점선 (:)
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.title('van der Pol Equation ($\mu=1$)')
plt.grid(True)
plt.show()
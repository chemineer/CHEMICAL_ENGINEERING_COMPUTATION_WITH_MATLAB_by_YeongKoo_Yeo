import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 1. 미분 방정식 시스템 정의
def tank_system(t, y):
    # y[0]=x1, y[1]=x2, y[2]=x3, y[3]=y4
    dy1_dt = (2.25 - 27 * y[0]) / y[3]
    dy2_dt = 15 * (y[0] - y[1]) / 20
    dy3_dt = 15 * (y[1] - y[2]) / 20
    dy4_dt = 12
    return [dy1_dt, dy2_dt, dy3_dt, dy4_dt]

# 2. 초기값 및 시간 범위 설정
tspan = (0, 2)
y0 = [0.15, 0.15, 0.15, 20]

# 3. ODE 해결 (해상도를 위해 t_eval 설정)
sol = solve_ivp(tank_system, tspan, y0, t_eval=np.linspace(0, 2, 200))

# 4. 결과 시각화
plt.figure(figsize=(8, 5))

# y1: 실선, y2: 점-선(.-), y3: 점선(:)
plt.plot(sol.t, sol.y[0], label='x_1')
plt.plot(sol.t, sol.y[1], '.-', label='x_2', markevery=10) # 마커가 너무 많지 않게 조절
plt.plot(sol.t, sol.y[2], ':', label='x_3')

plt.xlabel('t(min)')
plt.ylabel('x(t)')
plt.legend(loc='best')
plt.title('Well-Mixed Tanks System')
plt.grid(True, linestyle='--', alpha=0.7)

plt.show()
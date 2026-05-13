import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 1. ODE 정의: dy/dt = exp(-t)
def dy(t, y):
    return np.exp(-t)

# 2. 초기값 및 설정
tspan = (0, 1)  # 시간 범위
y0 = [-1]       # 초기값 (리스트 또는 배열 형태여야 함)

# 3. ODE 해결 (MATLAB의 ode45에 해당)
# solve_ivp는 기본 solver가 'RK45'입니다.
sol = solve_ivp(dy, tspan, y0, t_eval=np.linspace(0, 1, 100))

# 4. 결과 시각화
plt.plot(sol.t, sol.y[0])
plt.grid(True)
plt.xlabel('t')
plt.ylabel('y(t)')
plt.title('Solution of ODE: dy/dt = exp(-t)')
plt.show()
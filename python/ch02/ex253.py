import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 초기값 설정
t_span = (0, 20)
y0 = [0.05, 5]
k = 0.3
K = 1e-6

# 2. 강성 미분 방정식 정의 (Stiff ODE)
def stiff_ode(t, y):
    y1, y2 = y
    # 분모에 매우 작은 값 K가 있어 수치적으로 민감할 수 있음
    rate = k * y1 * y2 / (K + y2)
    dy1dt = rate
    dy2dt = -0.75 * rate
    return [dy1dt, dy2dt]

# 3. ODE 풀이 
# method='BDF' (Backward Differentiation Formula)는 MATLAB의 ode15s와 유사한 강성 솔버입니다.
sol = solve_ivp(stiff_ode, t_span, y0, method='BDF', t_eval=np.linspace(0, 20, 200))

# 4. 결과 시각화
plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0], label='B(t) (Biomass)', linewidth=1.5)
plt.plot(sol.t, sol.y[1], ':', label='S(t) (Substrate)', linewidth=1.5)

plt.grid(True)
plt.xlabel('t(min)')
plt.ylabel('y(t)')
plt.title('Solution of Stiff Differential Equations')
plt.legend(loc='best')
plt.show()
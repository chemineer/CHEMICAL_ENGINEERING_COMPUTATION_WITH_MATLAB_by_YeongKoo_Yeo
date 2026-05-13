import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 물성치 설정
Ta = 25          # 주변 온도 (deg.C)
L = 0.01         # 길이 (m)
Ac = 6.28e-5     # 대류 표면적 (m^2)
Ak = 3.14e-6     # 전도 단면적 (m^2)
hc = 20          # 대류 계수 (W/m^2.K)

rho = 2707       # 밀도 (kg/m^3)
k = 220          # 열전도도 (W/m.K)
Cp = 896         # 비열 (J/kg.K)

# 2. 계산된 매개변수
V = Ak * L       # 부피
Ch = rho * Cp * V   # 열용량
Rc = 1 / (hc * Ac)  # 대류 열저항
Rk = L / (k * Ak)   # 전도 열저항

T0 = Ta          # 초기 온도 (t=0일 때 핀의 온도)
Tb = 100         # 기부(Base) 온도 (deg.C)
tspan = (0, 8)   # 시간 범위

# 3. 미분 방정식 정의: dT/dt = f(t, T)
def dT_dt(t, T):
    term1 = (1 / (Ch * Rc)) * Ta
    term2 = (1 / (Ch * Rk)) * Tb
    term3 = (1 / (Ch * Rc) + 1 / (Ch * Rk)) * T
    return term1 + term2 - term3

# 4. ODE 풀이 (MATLAB의 ode45와 유사한 solve_ivp 사용)
sol = solve_ivp(dT_dt, tspan, [T0], method='RK45', t_eval=np.linspace(0, 8, 100))

# 5. 결과 시각화
plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0])
plt.xlabel('t(s)')
plt.ylabel('T(deg.C)')
plt.title('Pin Temperature Profile')
plt.grid(True)
plt.show()

# 6. 정상 상태(Steady-state) 온도 계산 및 출력
Ts = ((1 / (Ch * Rc)) * Ta + (1 / (Ch * Rk)) * Tb) / (1 / (Ch * Rc) + 1 / (Ch * Rk))
print(f"At steady-state, T(pin) = {Ts:.4f} deg.C")
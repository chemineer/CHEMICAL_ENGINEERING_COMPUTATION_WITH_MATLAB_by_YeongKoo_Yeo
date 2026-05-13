import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# 1. 데이터 및 물성치 설정
Ta = 25          # 주변 온도 (deg.C)
L0 = 0.01        # 전체 길이 (m)
D = 0.002         # 직경 (m)
L = L0 / 5       # 요소당 길이 (m)
Ac = np.pi * D * L  # 요소의 대류 표면적 (m^2)
Ak = 3.14e-6     # 전도 단면적 (m^2)
hc = 20          # 대류 계수 (W/m^2.K)

rho = 2707       # 밀도 (kg/m^3)
k = 220          # 열전도도 (W/m.K)
Cp = 896         # 비열 (J/kg.K)

# 2. 계산된 매개변수
V = Ak * L       # 부피
Ch = rho * Cp * V   # 열용량
Rc = 1 / (hc * Ac)  # 대류 열저항
Rce = 1 / (hc * Ak) # 끝부분 대류 열저항
Rk = L / (k * Ak)   # 요소 간 전도 열저항
Rk0 = L / (2 * k * Ak) # 기부(Base)와 첫 번째 요소 간 전도 열저항

T0 = np.full(5, Ta) # 초기 온도 (모두 Ta로 설정)
Tb = 100            # 기부 온도 (deg.C)
tspan = (0, 3)      # 시간 범위

# 3. 미분 방정식 시스템 정의 (dT/dt)
def system_dynamics(t, T):
    T1, T2, T3, T4, T5 = T
    
    dT1 = (-(1/Rk0 + 1/Rk + 1/Rc)*T1 + T2/Rk + Ta/Rc + Tb/Rk0) / Ch
    dT2 = (T1/Rk - (2/Rk + 1/Rc)*T2 + T3/Rk + Ta/Rc) / Ch
    dT3 = (T2/Rk - (2/Rk + 1/Rc)*T3 + T4/Rk + Ta/Rc) / Ch
    dT4 = (T3/Rk - (2/Rk + 1/Rc)*T4 + T5/Rk + Ta/Rc) / Ch
    dT5 = (T4/Rk - (1/Rk + 1/Rc + 1/Rce)*T5 + (1/Rc + 1/Rce)*Ta) / Ch
    
    return [dT1, dT2, dT3, dT4, dT5]

# 4. ODE 풀이 (MATLAB의 ode45 대응)
sol = solve_ivp(system_dynamics, tspan, T0, method='RK45', t_eval=np.linspace(0, 3, 300))

# 5. 결과 시각화
plt.figure(figsize=(10, 6))
styles = ['-', ':', '-.', '--', '.']
labels = ['T_1', 'T_2', 'T_3', 'T_4', 'T_5']

for i in range(5):
    plt.plot(sol.t, sol.y[i], styles[i], label=labels[i])

plt.xlabel('t(s)')
plt.ylabel('T(deg.C)')
plt.title('Temperature Profiles of Pin Elements')
plt.legend()
plt.grid(True)
plt.show()

# 6. 정상 상태(Steady-state) 온도 계산 (fsolve 사용)
def steady_state_eq(x):
    # 미분 방정식의 분자 부분이 0이 되는 지점을 찾음
    return [
        -(1/Rk0 + 1/Rk + 1/Rc)*x[0] + x[1]/Rk + Ta/Rc + Tb/Rk0,
        x[0]/Rk - (2/Rk + 1/Rc)*x[1] + x[2]/Rk + Ta/Rc,
        x[1]/Rk - (2/Rk + 1/Rc)*x[2] + x[3]/Rk + Ta/Rc,
        x[2]/Rk - (2/Rk + 1/Rc)*x[3] + x[4]/Rk + Ta/Rc,
        x[3]/Rk - (1/Rk + 1/Rc + 1/Rce)*x[4] + (1/Rc + 1/Rce)*Ta
    ]

x0 = np.full(5, Ta)
Ts = fsolve(steady_state_eq, x0)

print("Steady-state temperature:")
for i, temp in enumerate(Ts, 1):
    print(f"Pin element {i} = {temp:.4f} deg.C")
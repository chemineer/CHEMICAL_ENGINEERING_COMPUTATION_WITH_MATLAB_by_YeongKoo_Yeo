import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 파라미터 설정
rho = 980
V = 2.5
Cp = 1.6
w = 250
taut = 3.6
theta = 1.2
Tr = 75
Ti_initial = 50
Kc = 60
K = 1 / (Cp * w)
tau = rho * V / w
Qs = (Tr - Ti_initial) / K

# 2. 미분 방정식 정의 (매트랩의 hpf 함수 대응)
def hpf(t, z, tau, taut, K, Kc, theta, Qs, Tr):
    # z[0] = T (Tank), z[1] = T1 (Outlet), z[2] = Tt (Thermocouple)
    T, T1, Tt = z
    
    # 시간 t가 10분 이상일 때 유입 온도 Ti가 50에서 30으로 변화
    Ti = 50 if t < 10 else 30
    
    # 제어 신호(열량) Q 계산 (비례 제어)
    Q = Qs + Kc * (Tr - Tt)
    
    # 상태 방정식 계산
    dTdt = (Ti - T) / tau + K * Q / tau
    dT1dt = 2 * (T - T1 - (theta / 2) * dTdt) / theta
    dTtdt = (T1 - Tt) / taut
    
    return [dTdt, dT1dt, dTtdt]

# 3. 초기 조건 및 시뮬레이션 시간 설정
z0 = [Tr, Tr, Tr]  # 초기 온도 설정
t_span = (0, 60)   # 0분에서 60분까지
t_eval = np.linspace(0, 60, 1000) # 그래프를 위한 시간 간격

# 4. ODE 풀이 수행 (ode45와 유사한 RK45 방식)
sol = solve_ivp(
    fun=hpf,
    t_span=t_span,
    y0=z0,
    args=(tau, taut, K, Kc, theta, Qs, Tr),
    t_eval=t_eval,
    method='RK45'
)

# 5. 결과 시각화 (매트랩 plot 함수 대응)
plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='Tank (T)')
plt.plot(sol.t, sol.y[2], ':', label='Thermocouple (Tt)')
plt.plot(sol.t, sol.y[1], '.-', label='Outlet (T1)', markevery=40)

plt.xlabel('t(min)')
plt.ylabel('T(C)')
plt.title('CSTH P-Control Simulation (Example 9.20_1)')
plt.legend()
plt.grid(True)
plt.show()
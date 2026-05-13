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
tauI = 1.8
K = 1 / (Cp * w)
tau = rho * V / w
Qs = (Tr - Ti_initial) / K
Kc = 0 # 매트랩 코드의 Kc 값

# 2. 미분 방정식 정의 (매트랩의 htf 함수)
def htf(t, z, tau, taut, tauI, K, Kc, theta, Qs, Tr):
    # z[0]=T, z[1]=T1, z[2]=Tt, z[3]=Ce
    T, T1, Tt, Ce = z
    
    # 시간 조건에 따른 Ti 변화
    Ti = 50 if t < 10 else 30
    
    # 제어기 출력 Q 계산
    Q = Qs + (Kc / tauI) * Ce + Kc * (Tr - Tt)
    
    # 미분식 계산
    dTdt = (Ti - T) / tau + K * Q / tau
    dT1dt = 2 * (T - T1 - (theta / 2) * dTdt) / theta
    dTtdt = (T1 - Tt) / taut
    dCedt = Tr - Tt
    
    return [dTdt, dT1dt, dTtdt, dCedt]

# 3. 초기 조건 및 시뮬레이션 설정
z0 = [Tr, Tr, Tr, 0] # 초기값
t_span = (0, 60)    # tv = [0 60]
t_eval = np.linspace(0, 60, 600) # 그래프를 매끄럽게 그리기 위한 시간축

# 4. ODE 풀이 (solve_ivp 사용)
sol = solve_ivp(
    fun=htf, 
    t_span=t_span, 
    y0=z0, 
    args=(tau, taut, tauI, K, Kc, theta, Qs, Tr),
    t_eval=t_eval,
    method='RK45' # 매트랩의 ode45와 대응
)

# 5. 결과 시각화
plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='Tank (T)')
plt.plot(sol.t, sol.y[2], ':', label='Thermocouple (Tt)')
plt.plot(sol.t, sol.y[1], '.-', label='Outlet (T1)', markevery=20)

plt.xlabel('t(min)')
plt.ylabel('T(C)')
plt.title('CSTH Temperature Control Simulation')
plt.legend(loc='best')
plt.grid(True)
plt.show()
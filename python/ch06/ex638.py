import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
F = 1
Sf = 100

# 2. 미분 방정식 정의 (dx/dt)
def dx_dt(t, x):
    # x[0]=V (부피), x[1]=x (세포), x[2]=S (기질), x[3]=P (생성물)
    V, cell, S, P = x
    
    # 각 상태 변수의 변화율 계산
    # dV/dt
    dV = F
    # dx/dt (세포 성장)
    dx = (0.408 * S * np.exp(-0.028 * P) / (0.22 + S) - F / V) * cell
    # dS/dt (기질 소비)
    dS = -4.08 * S * cell * np.exp(-0.028 * P) / (0.22 + S) + F * (Sf - S) / V
    # dP/dt (생성물 형성)
    dP = S * cell * np.exp(-0.015 * P) / (0.44 + S) - F * P / V
    
    return [dV, dx, dS, dP]

# 3. 초기 조건 및 시간 구간 설정
tint = [0, 20]
x0 = [1, 0.2, 100, 0] # [V0, x0, S0, P0][cite: 17]

# 4. ODE 풀이 (ode45에 해당하는 RK45 방식 사용)[cite: 17]
sol = solve_ivp(dx_dt, tint, x0, method='RK45', dense_output=True, t_eval=np.linspace(0, 20, 200))

t = sol.t
V_val = sol.y[0]
x_val = sol.y[1]
S_val = sol.y[2]
P_val = sol.y[3]

# 5. 비성장 속도(mu) 및 비생산 속도(phi) 계산[cite: 17]
mu = 0.408 * S_val * np.exp(-0.028 * P_val) / (0.22 + S_val)
phi = S_val * np.exp(-0.015 * P_val) / (0.44 + S_val)

# 6. 결과 시각화[cite: 17]
plt.figure(figsize=(12, 5))

# Subplot 1: 농도 프로파일[cite: 17]
plt.subplot(1, 2, 1)
plt.plot(t, x_val, label='x')
plt.plot(t, S_val, '--', label='S')
plt.plot(t, P_val, ':', label='P')
plt.xlabel('t(hr)')
plt.ylabel('Concentration(g/liter)')
plt.grid(True)
plt.legend(loc='best')

# Subplot 2: 속도 프로파일 (mu 및 phi)[cite: 17]
plt.subplot(1, 2, 2)
plt.plot(t, mu, label='$\mu$')
plt.plot(t, phi, '--', label='$\pi$')
plt.xlabel('t(hr)')
plt.ylabel('$\mu$ and $\pi$')
plt.grid(True)
plt.legend(loc='best')

plt.tight_layout()
plt.show()
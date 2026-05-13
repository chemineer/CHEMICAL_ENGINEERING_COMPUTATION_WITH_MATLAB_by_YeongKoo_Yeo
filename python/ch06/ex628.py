import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
u = 8
Ap = 12
k = 30
Ka = 5
Ca0 = 0.2
kds = 17.5
kdp = 140

# 2. ODE 시스템 정의
def dydz(z, y):
    # y[0]=Xa1, y[1]=Xa2, y[2]=Xa3, y[3]=a2, y[4]=a3
    Xa1, Xa2, Xa3, a2, a3 = y
    
    # 각 미분 방정식 구현
    d_Xa1 = (1 / (1 + Ap * np.sqrt(z/u))) * k * Ca0 * (1 - Xa1) / (1 + Ka * Ca0 * (1 - Xa1)) / u
    d_Xa2 = a2 * k * Ca0 * (1 - Xa2) / (1 + Ka * Ca0 * (1 - Xa2)) / u
    d_Xa3 = a3 * k * Ca0 * (1 - Xa3) / (1 + Ka * Ca0 * (1 - Xa3)) / u
    d_a2 = -kds * (a2**2) / u
    d_a3 = -kdp * Ca0 * Xa3 * a3 / u
    
    return [d_Xa1, d_Xa2, d_Xa3, d_a2, d_a3]

# 3. 초기 조건 및 적분 구간 설정[cite: 5]
zspan = [0, 6]
y0 = [0, 0, 0, 1, 1]

# 4. ODE 풀이[cite: 5]
sol = solve_ivp(dydz, zspan, y0, method='RK45', t_eval=np.linspace(0, 6, 100))

z = sol.t
y = sol.y

# 5. 시각화[cite: 5]
plt.figure(figsize=(12, 5))

# Subplot 1: 전환율 (Xa1, Xa2, Xa3)
plt.subplot(1, 2, 1)
plt.plot(z, y[0], label='X_{A1}', linestyle='-')
plt.plot(z, y[1], label='X_{A2}', linestyle=':')
plt.plot(z, y[2], label='X_{A3}', linestyle='--')
plt.xlabel('z(m)')
plt.ylabel('X_A')
plt.legend()
plt.grid(True)

# Subplot 2: 촉매 활성도 (as, ap)
plt.subplot(1, 2, 2)
plt.plot(z, y[3], label='a_s', linestyle='-')
plt.plot(z, y[4], label='a_p', linestyle=':')
plt.xlabel('z(m)')
plt.ylabel('a')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
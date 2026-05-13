import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 데이터 설정
q = 0.12     # 유량
Ac = 0.26    # 단면적
k = 2.1      # 반응 속도 상수
kr = 0.98    # 흡착 평형 상수 관련 파라미터
Ca0 = 1      # 초기 농도

# (1) 농도 및 전환율 프로파일 (z에 따른 변화)
z_span = (0, 0.5)
z_eval = np.linspace(0, 0.5, 100) # 그래프를 매끄럽게 그리기 위한 구간 설정

# ODE 정의: dCa/dz
def dC_dz(z, Ca):
    return -(Ac / q) * k * Ca / np.sqrt(1 + kr * Ca**2)

# ODE 풀이
sol1 = solve_ivp(dC_dz, z_span, [Ca0], t_eval=z_eval)
z = sol1.t
Ca = sol1.y[0]
x = (Ca0 - Ca) / Ca0 # 전환율 계산

# 그래프 출력
plt.figure(figsize=(12, 5))

# 농도 그래프
plt.subplot(1, 2, 1)
plt.plot(z, Ca)
plt.grid(True)
plt.xlabel('z(m)')
plt.ylabel('Ca(mol/m^3)')
plt.title('Concentration Profile')

# 전환율 그래프
plt.subplot(1, 2, 2)
plt.plot(z, x)
plt.grid(True)
plt.xlabel('z(m)')
plt.ylabel('x(conversion)')
plt.title('Conversion Profile')

plt.tight_layout()
plt.show()

# (2) 반응기 부피 계산 (전환율 80% 달성 기준)
V0 = 0
xd = 0.8
x_span = (0, xd)

# ODE 정의: dV/dx
def dV_dx(x, V):
    return (q / k) * np.sqrt(1 + kr * Ca0**2 * (1 - x)**2) / (1 - x)

# ODE 풀이
sol2 = solve_ivp(dV_dx, x_span, [V0])
V_final = sol2.y[0][-1]

print(f"Reactor volume for 80 percent conversion = {V_final:.6f} m^3")
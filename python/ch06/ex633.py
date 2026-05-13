import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
Fa0 = 2.26e-5
R = 8.314
P0 = 1641
T0 = 698
E = 24000
Vf = 1e-5

# 2. 기초 계산
Ct0 = P0 / (R * T0)
# 기체 상수 R의 단위에 맞춰 에너지 항 계산 (MATLAB 코드의 1.987은 kcal/mol 단위 상수로 추정됨)
k = 0.29 * np.exp(E / 1.987 * (1 / 500 - 1 / T0))

# 3. 미분 방정식 정의
# x[0]: Fa, x[1]: Fb, x[2]: Fc
def dxdv(v, x):
    total_F = np.sum(x)
    term = k * (Ct0**2) * (x[0] / total_F)**2
    
    dfa_dv = -term
    dfb_dv = term
    dfc_dv = term / 2
    return [dfa_dv, dfb_dv, dfc_dv]

# 4. 초기 조건 및 적분 구간 설정[cite: 12]
x0 = [Fa0, 0, 0]
v_span = [0, Vf]

# 5. ODE 풀이 (MATLAB의 ode23s와 유사한 Radau 또는 BDF 방식 사용)[cite: 12]
sol = solve_ivp(dxdv, v_span, x0, method='Radau', dense_output=True)

V = sol.t
Fa = sol.y[0, :]
Fb = sol.y[1, :]
Fc = sol.y[2, :]

# 6. 결과 출력[cite: 12]
print(f'Final values: Fa = {Fa[-1]:g}, Fb = {Fb[-1]:g}, Fc = {Fc[-1]:g}')

# 7. 시각화[cite: 12]
plt.figure(figsize=(8, 5))
plt.plot(V, Fa, label='F$_A$', linestyle='-')
plt.plot(V, Fb, label='F$_B$', linestyle=':')
plt.plot(V, Fc, label='F$_C$', linestyle='--')

plt.xlabel('V(dm$^3$)')
plt.ylabel('Molar flow rate (mol/min)')
plt.legend()
plt.axis('tight')
plt.grid(True)
plt.show()
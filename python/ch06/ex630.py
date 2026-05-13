import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from C2H6fun import C2H6fun  # C2H6fun.py 임포트

# 1. 데이터 설정
# 1: C2H6, 2: C2H4, 3: C2H2, 4: C3H6; 5: C3H8, 6: C4H6, 7: CH4, 8: H2
P = 3
T = 1073.15
Ri = 0.08314
F10 = 1300
Fs = 0.4 * F10

# 2. 초기 조건 및 구간 설정
Vspan = [0, 20000]
F0 = np.zeros(8)
F0[0] = F10  # 에탄 초기 유량[cite: 7]

# 3. ODE 시스템 풀이 (ode45에 해당)[cite: 7]
sol = solve_ivp(
    C2H6fun, 
    Vspan, 
    F0, 
    args=(P, T, Fs), 
    method='RK45',
    dense_output=True
)

V = sol.t
F = sol.y.T  # (n_points, 8) 형태

# 4. 각 지점에서의 농도(C) 계산[cite: 7]
n = len(V)
C = np.zeros((n, 8))
for k in range(n):
    Ft = np.sum(F[k, :])
    C[k, :] = F[k, :] * P / (Ft * Ri * T)

C1 = C[:, 0]  # C2H6 농도[cite: 7]
C2 = C[:, 1]  # C2H4 농도[cite: 7]

# 5. 결과 시각화[cite: 7]
plt.figure(figsize=(8, 6))
plt.plot(V, C1, label='C$_2$H$_6$', linestyle='-')
plt.plot(V, C2, label='C$_2$H$_4$', linestyle='--')
plt.xlabel('V(liter)')
plt.ylabel('C(mol/l)')
plt.legend()
plt.grid(True)
plt.show()

# 6. 전환율 및 최종 유량 출력[cite: 7]
x = (F10 - F[-1, 0]) / F10  # 에탄 전환율 계산[cite: 7]
print(f'Fractional conversion of ethane = {x:g}')
print('Flow rate of each component at the end of the reactor volume:')
print(f'C2H6 = {F[-1, 0]:g} mol/s, C2H4 = {F[-1, 1]:g} mol/s')
print(f'C2H2 = {F[-1, 2]:g} mol/s, C3H6 = {F[-1, 3]:g} mol/s')
print(f'C3H8 = {F[-1, 4]:g} mol/s, C4H6 = {F[-1, 5]:g} mol/s')
print(f'CH4 = {F[-1, 6]:g} mol/s,  H2 = {F[-1, 7]:g} mol/s')
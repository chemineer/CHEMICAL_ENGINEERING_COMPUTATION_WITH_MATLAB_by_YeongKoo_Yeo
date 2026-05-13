import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from exbco import exbco  # exbco.py 임포트

# 1. 데이터 설정
Ca0 = 9.3; Fa0 = 14.67; T1 = 360; T2 = 333; k1 = 31.1; K2 = 3.03
E = 65700; R = 8.314; dH = -6900; Cp0 = 158.889; Cpc = 28; Ua = 5000; m = 500
hx = 'cn'

# 2. 이분법 설정
criT = 1e-3
Taf_target = 310
errT = 10
Vi = [0, 5]

if hx == 'cn':
    Ta01 = 305.0
    Ta02 = 320.0
    Ta0m = (Ta01 + Ta02) / 2
    
    while errT >= criT:
        # 각 가정된 Ta0에 대한 ODE 풀이
        # 초기값 z0 = [Ta0, X=0, T=305]
        def solve_at(Ta_guess):
            sol = solve_ivp(exbco, Vi, [Ta_guess, 0, 305], 
                            args=(Ca0, Fa0, Cp0, Cpc, Ua, m, T1, T2, k1, K2, E, R, dH, hx))
            return sol.y[0, -1], sol # 최종 위치에서의 Ta값 반환

        Ta_end1, _ = solve_at(Ta01)
        Ta_endm, sol_m = solve_at(Ta0m)
        
        # 이분법 업데이트[cite: 6]
        if (Ta_end1 - Taf_target) * (Ta_endm - Taf_target) < 0:
            Ta02 = Ta0m
        else:
            Ta01 = Ta0m
        
        Ta0m = (Ta01 + Ta02) / 2
        errT = abs(Ta_endm - Taf_target)
    
    zm_y = sol_m.y
    V = sol_m.t
else:
    z0 = [310, 0, 305]
    sol = solve_ivp(exbco, Vi, z0, args=(Ca0, Fa0, Cp0, Cpc, Ua, m, T1, T2, k1, K2, E, R, dH, hx))
    zm_y = sol.y
    V = sol.t

# 3. 추가 변수 계산[cite: 6]
Ta = zm_y[0, :]
X = zm_y[1, :]
T = zm_y[2, :]

k = k1 * np.exp(E * (1/T1 - 1/T) / R)
Kc = K2 * np.exp(dH * (1/T2 - 1/T) / R)
ra = -k * Ca0 * (1 - (1 + 1/Kc) * X)
Xe = Kc / (1 + Kc)

# 4. 시각화[cite: 6]
plt.figure(figsize=(10, 8))

plt.subplot(2, 2, 1)
plt.plot(V, T, label='T')
plt.plot(V, Ta, '--', label='T_a')
plt.xlabel('V'), plt.ylabel('T(K)'), plt.legend()

plt.subplot(2, 2, 2)
plt.plot(V, X, label='X')
plt.plot(V, Xe, '--', label='X_e')
plt.xlabel('V'), plt.ylabel('X, X_e'), plt.legend()

plt.subplot(2, 2, 3)
plt.plot(V, -ra)
plt.xlabel('V'), plt.ylabel('-r_A')

plt.tight_layout()
plt.show()

# 5. 결과 출력[cite: 6]
print(f"Conversion (X) and equilibrium conversion (Xe): Xf = {X[-1]:g}, Xef = {Xe[-1]:g}")
print(f"Final T and Ta: Tf = {T[-1]:g}, Taf = {Ta[-1]:g}")
print(f"Final reaction rate: raf = {-ra[-1]:g}")
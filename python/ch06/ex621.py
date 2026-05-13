import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 미분 방정식 정의 (단열 PFR에서의 전환율 변화)
def fad(V, X, Ca0, Fa0, T1, T2, k1, K2, E, R, dH):
    # X는 현재의 전환율 (Conversion)
    # 단열 조건에 따른 온도 계산 식
    T = 330 + 43.4265 * X
    
    # 아레니우스 식에 따른 반응 속도 상수 k와 평형 상수 Kc 계산
    k = k1 * np.exp(E * (1/T1 - 1/T) / R)
    Kc = K2 * np.exp(dH * (1/T2 - 1/T) / R)
    
    # 반응 속도 ra 및 전환율 변화율 dX/dV 계산
    ra = -k * Ca0 * (1 - (1 + 1/Kc) * X)
    dX_dV = -ra / Fa0
    
    return dX_dV

# 2. 데이터 및 파라미터 설정[cite: 19]
Ca0 = 9.3      # 초기 농도
Fa0 = 146.7    # 초기 몰 유량
T1 = 360       # k에 대한 기준 온도
T2 = 333       # K2에 대한 기준 온도
k1 = 31.1      # 기준 온도 T1에서의 속도 상수
K2 = 3.03      # 기준 온도 T2에서의 평형 상수
E = 65700      # 활성화 에너지 (J/mol)
R = 8.314      # 기체 상수
dH = -6900     # 반응열 (J/mol)
Vspan = (0, 5) # 반응기 부피 구간 (0 to 5)
X0 = [0]       # 초기 전환율[cite: 19]

# 3. ODE 풀기 (solve_ivp 사용)[cite: 19]
V_eval = np.linspace(0, 5, 100)
sol = solve_ivp(fad, Vspan, X0, t_eval=V_eval, args=(Ca0, Fa0, T1, T2, k1, K2, E, R, dH))

# 4. 결과 가공 (온도, 반응속도, 평형 전환율 계산)[cite: 19]
V = sol.t
X = sol.y[0]
T = 330 + 43.4265 * X
k = k1 * np.exp(E * (1/T1 - 1/T) / R)
Kc = K2 * np.exp(dH * (1/T2 - 1/T) / R)
ra = -k * Ca0 * (1 - (1 + 1/Kc) * X)
Xe = Kc / (1 + Kc) # 평형 전환율[cite: 19]

# 5. 결과 출력[cite: 19]
print(f"Conversion (X) and equilibrium conversion (Xe): Xf = {X[-1]:.4f}, Xef = {Xe[-1]:.4f}")
print(f"Final temperature: Tf = {T[-1]:.4f}")
print(f"Final reaction rate: raf = {-ra[-1]:.4f}")

# 6. 시각화[cite: 19]
plt.figure(figsize=(12, 10))

# 온도 그래프[cite: 19]
plt.subplot(2, 2, 1)
plt.plot(V, T)
plt.xlabel('V')
plt.ylabel('T(K)')
plt.grid(True)

# 반응 속도 그래프[cite: 19]
plt.subplot(2, 2, 2)
plt.plot(V, -ra)
plt.xlabel('V')
plt.ylabel('-r_A')
plt.grid(True)

# 전환율 및 평형 전환율 그래프[cite: 19]
plt.subplot(2, 2, 3)
plt.plot(V, X, label='X')
plt.plot(V, Xe, '--', label='Xe')
plt.xlabel('V')
plt.ylabel('X, X_e')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
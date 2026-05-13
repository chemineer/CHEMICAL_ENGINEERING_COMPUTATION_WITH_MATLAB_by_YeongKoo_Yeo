import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from batdat import *  # batdat.py에서 데이터 함수 임포트[cite: 13]

# 1. 데이터 로드
p = batdat()  # 파라미터 딕셔너리 획득[cite: 13]

# 2. 미분 방정식 시스템 정의
def batrxn(t, X):
    # X[0] = Ca, X[1] = Cb, X[2] = T(deg.C)
    Ca, Cb, T = X
    
    # 에너지 수지 관련 계수 계산
    a1 = p["dH1"] / (p["rho"] * p["Cp"])
    a2 = p["dH2"] / (p["rho"] * p["Cp"])
    a3 = p["Uj"] * p["AjV"] / (p["rho"] * p["Cp"])
    a4 = p["Uc"] * p["AcV"] / (p["rho"] * p["Cp"])
    
    # 속도 상수 k 계산 (Arrhenius 식)[cite: 12]
    # 온도는 섭씨를 켈빈으로 변환 (273.15 + T)[cite: 12]
    k1 = p["A1"] * np.exp(-p["E1"] / (p["R"] * (273.15 + T)))
    k2 = p["A2"] * np.exp(-p["E2"] / (p["R"] * (273.15 + T)))
    
    # dX/dt 정의[cite: 12]
    dCa = -k1 * Ca**2
    dCb = k1 * Ca**2 - k2 * Cb
    dT = a1 * k1 * Ca**2 + a2 * k2 * Cb + a3 * (p["Ts"] - T) - a4 * (T - p["Tc"])
    
    return [dCa, dCb, dT]

# 3. 초기 조건 및 시간 범위 설정[cite: 12]
X0 = [1.5, 0, 50]  # 초기 Ca=1.5, Cb=0, T=50도[cite: 12]
t_span = (0, 6000) # 0초부터 6000초까지[cite: 12]
t_eval = np.linspace(0, 6000, 1000) # 그래프 출력을 위한 시간 지점

# 4. ODE 풀기 (MATLAB의 ode45와 유사한 solve_ivp 사용)[cite: 12]
sol = solve_ivp(batrxn, t_span, X0, t_eval=t_eval)

# 5. 결과 시각화[cite: 12]
plt.figure(figsize=(12, 5))

# 농도 그래프 (Ca, Cb)[cite: 12]
plt.subplot(1, 2, 1)
plt.plot(sol.t, sol.y[0], label='$C_A(t)$')
plt.plot(sol.t, sol.y[1], '--', label='$C_B(t)$')
plt.xlabel('t(s)')
plt.ylabel('$C(kmol/m^3)$')
plt.legend()
plt.grid(True)

# 온도 그래프 (T)[cite: 12]
plt.subplot(1, 2, 2)
plt.plot(sol.t, sol.y[2], color='red')
plt.xlabel('t(s)')
plt.ylabel('T(deg.C)')
plt.grid(True)

plt.tight_layout()
plt.show()
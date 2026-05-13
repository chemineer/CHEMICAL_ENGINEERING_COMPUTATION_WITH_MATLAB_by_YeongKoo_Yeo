import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- 데이터 및 파라미터 설정 ---
v0 = 6          # 초기 공급 유량
v = 12          # 배출 유량
V = 200         # 반응기 부피
Ca0 = 2         # A의 초기 공급 농도[cite: 11]
Cb0 = 2         # B의 초기 공급 농도[cite: 11]
k = 0.5         # 반응 속도 상수[cite: 11]

# 시간 범위 설정: 0에서 20까지 0.01 간격[cite: 11]
t0, tf = 0, 20
tspan = (t0, tf)
# np.linspace는 시작과 끝을 정확히 지정하므로 t_span을 벗어나지 않습니다.
# 0부터 20까지 0.01 간격으로 하고 싶다면 약 2001개의 점이 필요합니다.
t_eval = np.linspace(t0, tf, 2001)

# 초기 조건: [Ca1, Cb1, Ca2, Cb2, Ca3, Cb3][cite: 11]
# 첫 번째 반응기에만 원료가 있고 나머지는 0인 상태[cite: 11]
C0 = [Ca0, Cb0, 0, 0, 0, 0]

# --- 미분 방정식 정의 (dC/dt) ---[cite: 11]
def dCdt(t, C):
    # C[0]:Ca1, C[1]:Cb1, C[2]:Ca2, C[3]:Cb2, C[4]:Ca3, C[5]:Cb3[cite: 11]
    
    # 반응기 1[cite: 11]
    dCa1_dt = (v0 * Ca0 - v * C[0] - k * V * C[0] * C[1]) / V
    dCb1_dt = (v0 * Cb0 - v * C[1] - k * V * C[0] * C[1]) / V
    
    # 반응기 2[cite: 11]
    dCa2_dt = (v * C[0] - v * C[2] - k * V * C[2] * C[3]) / V
    dCb2_dt = (v * C[1] - v * C[3] - k * V * C[2] * C[3]) / V
    
    # 반응기 3[cite: 11]
    dCa3_dt = (v * C[2] - v * C[4] - k * V * C[4] * C[5]) / V
    dCb3_dt = (v * C[3] - v * C[5] - k * V * C[4] * C[5]) / V
    
    return [dCa1_dt, dCb1_dt, dCa2_dt, dCb2_dt, dCa3_dt, dCb3_dt]

# --- 수치 적분 수행 (solve_ivp) ---
# MATLAB의 ode45와 동일한 RK45 알고리즘 사용[cite: 11]
sol = solve_ivp(dCdt, tspan, C0, t_eval=t_eval, method='RK45')

t = sol.t
# 결과 추출: C[0], C[2], C[4]가 각각 각 반응기의 A 농도[cite: 11]
Ca1 = sol.y[0]
Ca2 = sol.y[2]
Ca3 = sol.y[4]

# --- 시각화 ---[cite: 11]
plt.figure(figsize=(8, 6))
plt.plot(t, Ca1, '-', label='C_{A1}')   # 실선[cite: 11]
plt.plot(t, Ca2, ':', label='C_{A2}')   # 점선[cite: 11]
plt.plot(t, Ca3, '--', label='C_{A3}')  # 파선[cite: 11]

plt.xlabel('Time(min)')
plt.ylabel('Concentration(gmol/dm^3)')
plt.legend()
plt.title('Concentration Profiles in 3 Series CSTRs')
plt.grid(True)
plt.show()
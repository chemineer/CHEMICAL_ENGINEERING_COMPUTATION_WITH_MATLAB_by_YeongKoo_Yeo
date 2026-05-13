import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- 데이터 및 파라미터 설정 ---
v = 25          # 유량 (Flow rate)
k1 = 0.2        # 반응 속도 상수 1
k2 = 0.1        # 반응 속도 상수 2
Ca0 = 2.5       # A의 초기 농도
Cb0 = 0         # B의 초기 농도[cite: 10]
Cc0 = 0         # C의 초기 농도[cite: 10]

tspan = (0, 15)      # 시간 범위 (시작, 끝)[cite: 10]
t_eval = np.arange(0, 15.01, 0.01)  # MATLAB의 0:0.01:15와 대응[cite: 10]
C0 = [Ca0, Cb0, Cc0] # 초기 조건 배열[cite: 10]

# --- 미분 방정식 정의 (dC/dt) ---[cite: 10]
# C[0]=Ca, C[1]=Cb, C[2]=Cc
def dCdt(t, C):
    dAdt = -k1 * C[0]                 #[cite: 10]
    dBdt = k1 * C[0] - k2 * C[1]      #[cite: 10]
    dC_dt = k2 * C[1]                 #[cite: 10]
    return [dAdt, dBdt, dC_dt]

# --- 수치 적분 수행 (solve_ivp 사용) ---
# MATLAB의 ode45와 유사한 RK45 방식을 사용합니다.
sol = solve_ivp(dCdt, tspan, C0, t_eval=t_eval, method='RK45')

t = sol.t
Ca, Cb, Cc = sol.y

# --- 결과 분석 (최대 농도 및 시간) ---[cite: 10]
Cmax = np.max(Cb)                             # B의 최대 농도[cite: 10]
max_idx = np.argmax(Cb)                       # 최대 농도일 때의 인덱스
tmax = t[max_idx]                             # 최대 농도 도달 시간[cite: 10]
Vol = v * tmax                                # 필요한 반응기 부피[cite: 10]

print(f"B의 최대 농도 (Cmax): {Cmax:.6f}")
print(f"최대 농도 도달 시간 (tmax): {tmax:.6f}")
print(f"반응기 부피 (Vol): {Vol:.6f}")

# --- 시각화 ---[cite: 10]
plt.figure(figsize=(8, 6))
plt.plot(t, Ca, '-', label='A')   # 실선[cite: 10]
plt.plot(t, Cb, ':', label='B')   # 점선[cite: 10]
plt.plot(t, Cc, '--', label='C')  # 파선[cite: 10]

plt.xlabel('Time(min)')             #[cite: 10]
plt.ylabel('Concentration(mol/liter)') #[cite: 10]
plt.legend()                        #[cite: 10]
plt.title('Concentration Profiles in Series Reaction')
plt.grid(True)
plt.show()
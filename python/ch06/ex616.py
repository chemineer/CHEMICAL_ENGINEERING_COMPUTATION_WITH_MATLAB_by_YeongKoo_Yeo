import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 미분 방정식 시스템 정의 (Semi-batch Reactor)
def semibrx(t, x, k, v0, V0, Cb0):
    # x[0]=Ca, x[1]=Cb, x[2]=Cc, x[3]=Cd
    Ca, Cb, Cc, Cd = x
    
    # 반응 속도 rA (A + B -> C + D 가정)
    ra_val = -k * Ca * Cb
    
    # 현재 시간 t에서의 부피 V
    V = V0 + v0 * t
    
    # 각 성분에 대한 물질 수지 (Mole Balance)
    # dCa/dt = rA - (v0/V)*Ca
    # dCb/dt = rA + (Cb0 - Cb)*(v0/V)
    # dCc/dt = -rA - (v0/V)*Cc
    # dCd/dt = -rA - (v0/V)*Cd
    dCa = ra_val - (v0 * Ca / V)
    dCb = ra_val + (Cb0 - Cb) * v0 / V
    dCc = -ra_val - (v0 * Cc / V)
    dCd = -ra_val - (v0 * Cd / V)
    
    return [dCa, dCb, dCc, dCd]

# 2. 데이터 및 파라미터 설정[cite: 14]
k = 2.2      # 반응 속도 상수
v0 = 0.05    # 주입 유량
V0 = 5       # 초기 부피
Cb0 = 0.025  # 주입되는 B의 농도
Ca0 = 0.05   # 초기 A의 농도

# 3. 초기 조건 및 시간 범위 설정[cite: 14]
x0 = [Ca0, 0, 0, 0]  # 초기 농도: Ca=0.05, 나머지는 0
t_span = (0, 500)    # 0초부터 500초까지
t_eval = np.linspace(0, 500, 1000) # 그래프 출력을 위한 시간 지점

# 4. ODE 풀기 (solve_ivp 사용)[cite: 14]
# args를 통해 추가 파라미터를 전달합니다.
sol = solve_ivp(semibrx, t_span, x0, t_eval=t_eval, args=(k, v0, V0, Cb0))

# 5. 결과 가공 (부피 및 반응 속도 계산)[cite: 14]
t = sol.t
Ca, Cb, Cc, Cd = sol.y
V = V0 + v0 * t
rA = k * Ca * Cb  # 시간에 따른 반응 속도 계산[cite: 14]

# 6. 시각화[cite: 14]
plt.figure(figsize=(12, 5))

# 농도 변화 그래프[cite: 14]
plt.subplot(1, 2, 1)
plt.plot(t, Ca, label='$C_A$')
plt.plot(t, Cb, ':', label='$C_B$')
plt.plot(t, Cc, '.-', label='$C_C$', markersize=3)
plt.plot(t, Cd, '--', label='$C_D$')
plt.xlabel('t(sec)')
plt.ylabel('Concentration(mol/dm^3)')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)

# 반응 속도 변화 그래프[cite: 14]
plt.subplot(1, 2, 2)
plt.plot(t, rA, color='red')
plt.xlabel('t(sec)')
plt.ylabel('Reaction rate(mol/dm^3sec)')
plt.title('Reaction Rate Over Time')
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.show()
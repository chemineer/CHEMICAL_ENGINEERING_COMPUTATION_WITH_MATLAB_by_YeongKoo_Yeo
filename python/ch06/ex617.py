import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
v0 = 0.5      # 부피 유량 (Volumetric flow rate)
k = 0.3       # 반응 속도 상수
Vf = 2.4      # 최종 반응기 부피
V_span = (0, Vf)  # 부피 구간 (0 ~ 2.4)
C0 = [2, 0, 2]    # 초기 농도: Ca=2, Cb=0, Cc=2

# 2. 미분 방정식 정의 (부피 V에 따른 농도 C의 변화)
def dC_dV(V, C):
    # C[0]=Ca, C[1]=Cb, C[2]=Cc
    Ca = C[0]
    
    # 성분 A에 대한 물질 수지: dCa/dV = -ra / v0
    # 여기서 ra = 2 * k * Ca^2 (2차 반응 가정)
    dCa = -2 * k * (Ca**2) / v0
    
    # 성분 B에 대한 물질 수지: dCb/dV = k * Ca^2 / v0
    dCb = k * (Ca**2) / v0
    
    # 성분 C는 비활성 성분이므로 변화율 0
    dCc = 0
    
    return [dCa, dCb, dCc]

# 3. ODE 풀기 (solve_ivp 사용)[cite: 15]
# t_eval을 지정하여 매끄러운 그래프를 위한 데이터를 생성합니다.
V_eval = np.linspace(0, Vf, 100)
sol = solve_ivp(dC_dV, V_span, C0, t_eval=V_eval)

# 4. 결과 시각화[cite: 15]
plt.figure(figsize=(8, 6))
plt.plot(sol.t, sol.y[0], label='$C_A$')          # 성분 A 농도
plt.plot(sol.t, sol.y[1], '.-', label='$C_B$')    # 성분 B 농도 (마커 포함)
plt.plot(sol.t, sol.y[2], '--', label='$C_C$')    # 성분 C 농도 (파선)

plt.xlabel('V(length, m)')
plt.ylabel('Concentration(kmol/m^3)')
plt.title('Isothermal PFR Concentration Profile')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
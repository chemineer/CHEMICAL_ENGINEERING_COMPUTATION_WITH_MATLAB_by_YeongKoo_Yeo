import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 미분 방정식 시스템 정의 (PFR 모델)
def pfr_complex(V, F, k, Ct0):
    # F[0]=A, F[1]=B, F[2]=C, F[3]=D, F[4]=E, F[5]=F
    Ft = np.sum(F)  # 총 몰 유량 계산
    
    # 각 반응의 속도 정의 (MATLAB 소스 기반)
    # Ct0는 기준 농도, F(i)/Ft는 몰 분율을 의미함
    r1A = -k[0] * (Ct0**3) * F[0] * (F[1]**2) / (Ft**3)
    r2A = -k[1] * (Ct0**2) * F[0] * F[1] / (Ft**2)
    r3B = -k[2] * (Ct0**3) * F[1] * (F[2]**2) / (Ft**3)
    r4C = -k[3] * (Ct0**(5/3)) * F[2] * (F[0]**(2/3)) / (Ft**(5/3))
    
    # 각 성분에 대한 물질 수지 (dFi/dV)
    dFa = r1A + r2A + (2/3) * r4C
    dFb = (5/4) * r1A + (3/4) * r2A + r3B
    dFc = -r1A + 2 * r3B + r4C
    dFd = -(3/2) * r1A - (3/2) * r2A - r4C
    dFe = -(1/2) * r2A - (5/6) * r4C
    dFf = -2 * r3B
    
    return [dFa, dFb, dFc, dFd, dFe, dFf]

# 2. 데이터 및 초기 파라미터 설정[cite: 17]
k = [5, 2, 10, 5]               # 반응 속도 상수 k1, k2, k3, k4[cite: 17]
F0 = [10, 10, 0, 0, 0, 0]       # 초기 몰 유량 [Fa0, Fb0, Fc0, Fd0, Fe0, Ff0][cite: 17]
Vf = 10                         # 반응기 최종 부피 (liter)[cite: 17]
Ct0 = 2                         # 초기 총 농도[cite: 17]

# 3. ODE 풀기 (solve_ivp 사용)[cite: 17]
V_span = (0, Vf)
V_eval = np.linspace(0, Vf, 200) # 그래프 출력을 위한 부피 지점
sol = solve_ivp(pfr_complex, V_span, F0, t_eval=V_eval, args=(k, Ct0))

# 4. 결과 시각화[cite: 17]
plt.figure(figsize=(10, 7))
labels = ['F_A', 'F_B', 'F_C', 'F_D', 'F_E', 'F_F']
styles = ['--', ':', '-', '.-', '-', '-'] # MATLAB 스타일 반영[cite: 17]

for i in range(6):
    plt.plot(sol.t, sol.y[i], styles[i], label=labels[i])

plt.xlabel('V(liter)')
plt.ylabel('F_i(mol/min)')
plt.title('Molar Flow Rates in PFR with Complex Reactions')
plt.legend(loc='best')
plt.grid(True, linestyle='--', alpha=0.6)
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
Fa0 = 10
k = 0.7
kc = 0.2
Kc = 0.05
R = 8.314
P = 830.6
T = 500
Ct0 = P / (R * T) # 초기 총 농도
Vf = 500 # 최종 반응 부피
Vc = 400 # 분석하고자 하는 특정 부피

# 2. 미분 방정식 정의 (MATLAB의 mrf 함수)[cite: 13]
def mrf(v, x, k, kc, Kc, Ct0):
    # x[0] = Fa, x[1] = Fb, x[2] = Fc[cite: 13]
    Ft = np.sum(x) # 총 몰 유량
    
    # 반응 속도 rA 계산[cite: 13]
    rA = -k * Ct0 * (x[0] / Ft - Ct0 * x[1] * x[2] / (Kc * Ft**2))
    
    # 각 성분에 대한 dF/dV[cite: 13]
    dfa_dv = rA
    dfb_dv = -rA - kc * Ct0 * x[1] / Ft
    dfc_dv = -rA
    
    return [dfa_dv, dfb_dv, dfc_dv]

# 3. ODE 풀이 (ode45에 해당)[cite: 13]
x0 = [Fa0, 0, 0] # 초기 조건[cite: 13]
sol = solve_ivp(
    mrf, 
    [0, Vf], 
    x0, 
    args=(k, kc, Kc, Ct0), 
    method='RK45', 
    dense_output=True
)

V = sol.t
Fa, Fb, Fc = sol.y

# 4. 결과 출력 및 시각화[cite: 13]
print(f'Final values: Fa = {Fa[-1]:g}, Fb = {Fb[-1]:g}, Fc = {Fc[-1]:g}')

plt.figure(figsize=(8, 5))
plt.plot(V, Fa, label='F$_A$', linestyle='-')
plt.plot(V, Fb, label='F$_B$', linestyle=':')
plt.plot(V, Fc, label='F$_C$', linestyle='--')
plt.xlabel('V(dm$^3$)')
plt.ylabel('Molar flow rate (mol/min)')
plt.legend()
plt.grid(True)
plt.autoscale(enable=True, axis='x', tight=True)
plt.show()

# 5. 특정 부피 Vc(400)에서의 Fa 및 전환율 계산[cite: 13]
# solve_ivp의 dense_output 기능을 사용하여 더 정확한 값을 추출할 수 있습니다.
Fc_at_Vc = sol.sol(Vc)[0] 
Xa = (Fa0 - Fc_at_Vc) / Fa0 # 전환율 계산[cite: 13]

print(f'At V = {Vc:g}, Fa = {Fc_at_Vc:g} and the conversion of A = {Xa:g}')
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
V = 1        # 반응기 부피
F = 25       # 유량
Caf = 10     # 유입 농도
k1 = 50      # 반응 1 속도 상수 (A -> B)
k2 = 100     # 반응 2 속도 상수 (B -> D)
k3 = 10      # 반응 3 속도 상수 (2A -> C)

# 2. 미분 방정식 시스템 정의
def vandrxn(t, C):
    # C[0] = Ca, C[1] = Cb
    Ca, Cb = C
    
    # dCa/dt: A의 농도 변화율
    # -k1*Ca (B로 변함), -k3*Ca^2 (C로 변함), F/V*(Caf - Ca) (유입 및 유출)
    dCa = -k1 * Ca - k3 * (Ca**2) + F * (Caf - Ca) / V
    
    # dCb/dt: B의 농도 변화율
    # k1*Ca (A로부터 생성), -k2*Cb (D로 변함), -F/V*Cb (유출)[cite: 10]
    dCb = k1 * Ca - k2 * Cb - F * Cb / V
    
    return [dCa, dCb]

# 3. 초기 조건 및 시간 범위 설정[cite: 10]
C0 = [10, 0]      # 초기 농도: Ca=10, Cb=0
t_span = (0, 0.06) # 시간 범위: 0부터 0.06까지
t_eval = np.linspace(0, 0.06, 200) # 그래프를 그리기 위한 시간 지점들

# 4. ODE 풀기 (solve_ivp 사용)[cite: 10]
sol = solve_ivp(vandrxn, t_span, C0, t_eval=t_eval)

# 5. 결과 시각화[cite: 10]
plt.figure(figsize=(8, 6))
plt.plot(sol.t, sol.y[0], label='$C_A$', color='blue')
plt.plot(sol.t, sol.y[1], ':', label='$C_B$', color='red')

plt.xlabel('t')
plt.ylabel('C(t)')
plt.title('Van de Vusse Reaction in CSTR')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
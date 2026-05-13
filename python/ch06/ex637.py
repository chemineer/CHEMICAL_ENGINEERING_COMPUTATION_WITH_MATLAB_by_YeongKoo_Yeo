import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
D = 0.1      # 희석률 (Dilution rate)
Sf = 5       # 유입 기질 농도
y1 = 0.8     # 세포 성장 수율[cite: 16]
y2 = 0.7     # 생성물 형성 수율[cite: 16]
mum = 0.6    # 최대 비성장 속도[cite: 16]
Km = 0.28    # 반포화 상수[cite: 16]

# 2. 미분 방정식 정의[cite: 16]
def dzdt(t, z):
    # z[0] = x (세포 농도), z[1] = S (기질 농도), z[2] = P (생성물 농도)[cite: 16]
    x, S, P = z
    
    # 비성장 속도 (Monod 식 적용)
    mu = mum * S / (Km + S)
    
    # 각 성분의 시간 따른 농도 변화율[cite: 16]
    dxdt = -D * x + y1 * mu * x
    dSdt = -mu * x + D * (Sf - S)
    dPdt = y2 * mu * x - D * P
    
    return [dxdt, dSdt, dPdt]

# 3. 초기 조건 및 시간 구간 설정[cite: 16]
X0 = 0.03
S0 = 5
P0 = 0
z0 = [X0, S0, P0]
tint = [0, 30]

# 4. ODE 풀이 (MATLAB의 ode15s에 해당하는 BDF 방식 사용)[cite: 16]
sol = solve_ivp(
    dzdt, 
    tint, 
    z0, 
    method='BDF', 
    dense_output=True,
    t_eval=np.linspace(tint[0], tint[1], 300)
)

t = sol.t
x_val = sol.y[0]
S_val = sol.y[1]
P_val = sol.y[2]

# 5. 결과 시각화[cite: 16]
plt.figure(figsize=(8, 6))
plt.plot(t, x_val, label='x', linestyle='-')
plt.plot(t, S_val, label='S', linestyle='--')
plt.plot(t, P_val, label='P', linestyle=':')

plt.xlabel('t(hr)')
plt.ylabel('Concentration(g/ml)')
plt.legend(loc='best')
plt.grid(True)
plt.show()
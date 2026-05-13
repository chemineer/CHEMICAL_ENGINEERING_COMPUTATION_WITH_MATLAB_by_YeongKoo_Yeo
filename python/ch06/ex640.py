import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 구조 설정
class PData:
    def __init__(self):
        self.Fv = 0.0625
        self.S0 = 0.5
        self.m1 = 0.25
        self.m2 = 0.24
        self.K1 = 5e-4
        self.K2 = 4e8
        self.c1 = 3.3e-10
        self.c2 = 1.4e3

pdat = PData()

# 2. 미분 방정식 정의
def pred(t, x, pdat):
    S, N1, N2 = x
    
    S0, Fv, c1, c2 = pdat.S0, pdat.Fv, pdat.c1, pdat.c2
    m1, m2, K1, K2 = pdat.m1, pdat.m2, pdat.K1, pdat.K2
    
    # 미분 방정식 정의
    # dS/dt
    dS = Fv * (S0 - S) - (c1 * m1 * N1 * S) / (K1 + S)
    # dN1/dt
    dN1 = -Fv * N1 + (m1 * N1 * S) / (K1 + S) - (c2 * m2 * N1 * N2) / (K2 + N1)
    # dN2/dt
    dN2 = -Fv * N2 + (m2 * N1 * N2) / (K2 + N1)
    
    return [dS, dN1, dN2]

# 3. 초기 조건 및 시간 구간 설정
N10 = 1.3e9
N20 = 4e5
x0 = [pdat.S0, N10, N20]
tspan = [0, 1000]

# 4. ODE 풀이 (수치적 안정성을 위해 오차 한계 설정 추가)
sol = solve_ivp(
    pred, 
    tspan, 
    x0, 
    args=(pdat,), 
    method='BDF', 
    rtol=1e-6,  # 상대 오차 한계 설정
    atol=1e-9,  # 절대 오차 한계 설정
    t_eval=np.linspace(tspan[0], tspan[1], 1000)
)

t = sol.t
N1_val = sol.y[1]
N2_val = sol.y[2]

# 5. 결과 시각화
plt.figure(figsize=(10, 6))
plt.plot(t, np.log10(np.maximum(N1_val, 1e-1)), label='Bacteria(N$_1$)', linestyle='-')
plt.plot(t, np.log10(np.maximum(N2_val, 1e-1)), label='Amoeba(N$_2$)', linestyle='--')

plt.xlabel('t(hr)')
plt.ylabel('log$_{10}$(N$_1$) and log$_{10}$(N$_2$)')
plt.legend(loc='best')
plt.grid(True)
plt.title('Corrected Predator-Prey Model (E.Coli - Amoeba)')
plt.show()
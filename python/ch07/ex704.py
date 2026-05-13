import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 초기 조건 및 파라미터 설정
x10, x20, x30 = 0.5, 0.4, 0.35
x40, x50 = 0.0, 0.0
V = 1200
tau = 240
Sa = 1
Sd = 0.4
rho = 1414.7

tspan = [0, 150]
x0 = [x10, x20, x30, x40, x50]

# 2. 미분 방정식 정의 (ruf 함수)
def ruf(t, x, V, rho, tau, Sa, Sd):
    # x[0,1,2] = D (입자 직경), x[3] = Cas, x[4] = Cds
    D = x[:3]
    Cas = x[3]
    Cds = x[4]
    
    fd = np.zeros(5)
    sum0 = 0.0
    sum1 = 0.0
    
    for i in range(3):
        kL = 1.2 / D[i] # 물질 전달 계수
        Sw = 0
        
        # 입자 크기에 따른 용해 속도 및 상태 결정[cite: 13]
        if D[i] > 0.3:
            fd[i] = -2 * kL * (Sa - Cas) / rho
            Sw = 1
        elif 0.3 >= D[i] >= 1e-5:
            fd[i] = -2 * kL * (Sd - Cds) / rho
        else:
            fd[i] = 0
            
        sum0 += Sw * kL * D[i]
        sum1 += (1 - Sw) * kL * D[i]
    
    # 농도 변화율 계산[cite: 13]
    fd[3] = np.pi * (Sa - Cas) * sum0 / V - Cas / tau
    fd[4] = np.pi * (Sd - Cds) * sum1 / V - Cas / tau # 원문 코드의 x(4)/tau 유지[cite: 13]
    
    return fd

# 3. ODE 풀이[cite: 13]
sol = solve_ivp(
    ruf, 
    tspan, 
    x0, 
    args=(V, rho, tau, Sa, Sd), 
    method='RK45', 
    t_eval=np.linspace(tspan[0], tspan[1], 500)
)

t = sol.t
D1, D2, D3 = sol.y[0], sol.y[1], sol.y[2]
Cas, Cds = sol.y[3], sol.y[4]

# 4. 결과 시각화[cite: 13]
plt.figure(figsize=(12, 5))

# Subplot 1: 입자 직경 변화[cite: 13]
plt.subplot(1, 2, 1)
plt.plot(t, D1, label='D$_1$')
plt.plot(t, D2, ':', label='D$_2$')
plt.plot(t, D3, '-.', label='D$_3$')
plt.xlabel('t(min)')
plt.ylabel('D(cm)')
plt.legend()
plt.grid(True)
plt.autoscale(enable=True, axis='x', tight=True)

# Subplot 2: 약물 농도 변화[cite: 13]
plt.subplot(1, 2, 2)
plt.plot(t, Cas, label='C$_{AS}$')
plt.plot(t, Cds, ':', label='C$_{DS}$')
plt.xlabel('t(min)')
plt.ylabel('C(mg/cm$^3$)')
plt.legend(loc='best')
plt.grid(True)
plt.autoscale(enable=True, axis='x', tight=True)

plt.tight_layout()
plt.show()
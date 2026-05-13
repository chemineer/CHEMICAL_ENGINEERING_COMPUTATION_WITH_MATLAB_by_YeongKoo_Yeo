import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from adfun import adfun  # adfun.py 모듈 임포트

# 1. 데이터 및 파라미터 설정
# P: 1.6 atm부터 5 atm까지 0.2 간격으로 설정 후 kPa로 변환
P_atm = np.arange(1.6, 5.1, 0.2)
P_kPa = P_atm * 101.325
Vspan = [0, 4]

# 초기 공급 유량 및 질소 유량 설정
FA0 = np.array([10, 20, 30, 35, 38.3])
FN2 = 38.3 - FA0

nP = len(P_kPa)
nF = len(FA0)

# 결과 저장을 위한 행렬 초기화
xc = np.zeros((nF, nP))
T_final = np.zeros((nF, nP))

# 2. 이중 루프를 이용한 시뮬레이션 수행[cite: 10]
for i in range(nF):
    for j in range(nP):
        # 초기 조건: [FA, FB, FC, T] (초기 온도 1035K)[cite: 10]
        X0 = [FA0[i], 0, 0, 1035]
        pf = [P_kPa[j], FN2[i]]
        
        # ODE 풀이[cite: 10]
        sol = solve_ivp(adfun, Vspan, X0, args=(pf,), method='RK45')
        
        # 최종 전환율 및 온도 저장[cite: 10]
        # xc = (초기 유량 - 최종 유량) / 초기 유량[cite: 10]
        xc[i, j] = (X0[0] - sol.y[0, -1]) / X0[0]
        T_final[i, j] = sol.y[3, -1]

# 3. 결과 시각화[cite: 10]
plt.figure(figsize=(14, 6))
markers = ['o', '*', 'x', 'd', 'v']
labels = ['F_A0=10', 'F_A0=20', 'F_A0=30', 'F_A0=35', 'F_A0=38.3']

# Subplot 1: 압력에 따른 최종 온도[cite: 10]
plt.subplot(1, 2, 1)
for i in range(nF):
    plt.plot(P_atm, T_final[i, :], marker=markers[i], linestyle='None', label=labels[i])
plt.xlabel('P(atm)')
plt.ylabel('T(K)')
plt.legend()
plt.grid(True)

# Subplot 2: 압력에 따른 최종 전환율[cite: 10]
plt.subplot(1, 2, 2)
for i in range(nF):
    plt.plot(P_atm, xc[i, :], marker=markers[i], linestyle='None', label=labels[i])
plt.xlabel('P(atm)')
plt.ylabel('x_A')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
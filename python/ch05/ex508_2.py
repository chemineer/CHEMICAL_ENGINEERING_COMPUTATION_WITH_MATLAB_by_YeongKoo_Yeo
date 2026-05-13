import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from functools import partial
from vwfun import vwfun

# 1. 초기 파라미터 정의
T = 60
rf = 0.00015
dz = 300
dP = -150
D_list = np.array([4.026, 5.047, 6.065, 7.981])
L_list = np.arange(500, 10500, 500) # 500부터 10000까지 500 단위

# 결과 저장을 위한 배열 (행: L, 열: D)
v_results = np.zeros((len(L_list), len(D_list)))
q_results = np.zeros((len(L_list), len(D_list)))

# 2. 이중 루프 계산
for i, D in enumerate(D_list):
    for j, L in enumerate(L_list):
        # 고정된 인자를 partial로 전달
        func = partial(vwfun, T=T, L=L, D=D, rf=rf, dz=dz, dP=dP)
        
        # 해 찾기 (초기값 10)
        v = fsolve(func, 10)[0]
        
        # 결과 저장
        v_results[j, i] = v
        q_results[j, i] = (7.481 * 60) * (np.pi * v * (D / 12)**2) / 4

# 3. 그래프 그리기
markers = ['o', '*', '+', 'd']
labels = ['D=4in', 'D=5in', 'D=6in', 'D=8in']

# Figure 1: Velocity
plt.figure(figsize=(8, 5))
for i in range(len(D_list)):
    plt.plot(L_list, v_results[:, i], markers[i], label=labels[i])
plt.legend(); plt.axis([500, 10000, 0, 20]); plt.xlabel('L (ft)'); plt.ylabel('Velocity, v(ft/s)')
plt.grid(True)

# Figure 2: Flow Rate
plt.figure(figsize=(8, 5))
for i in range(len(D_list)):
    plt.plot(L_list, q_results[:, i], markers[i], label=labels[i])
plt.legend(); plt.xlabel('L(ft)'); plt.ylabel('Flow rate, q(gpm)')
plt.grid(True)

plt.show()
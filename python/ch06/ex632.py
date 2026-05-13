import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pmmar import pmmar  # pmmar.py 임포트

# 1. 데이터 구조 설정 (pm 객체 대신 딕셔너리 사용)
pm = {
    'M0': 1.5e4,
    'MWm': 0.10013,
    'MWi': 0.077,
    'Mjp': 0.18781,
    'rhop': 1200,
    'Vms': 8.22e-4,
    'Vps': 7.7e-4,
    'Vis': 8.25e-4,
    'T': 350
}

# 2. 초기 조건 및 시간 구간 설정
tspan = [0, 6000]
M0 = pm['M0']
x0 = np.zeros(10)
x0[1] = M0  # MATLAB의 x0(2)는 인덱스 1

# 3. ODE 풀이 (MATLAB의 ode15s에 해당하는 BDF 방식 사용)
sol = solve_ivp(
    pmmar, 
    tspan, 
    x0, 
    args=(pm,), 
    method='BDF', 
    dense_output=True
)

t = sol.t
x = sol.y.T  # (n_points, 10)

# 4. 결과 추출 및 변환[cite: 11]
I, M, Q = x[:, 0], x[:, 1], x[:, 9]
L0, L1, L2 = x[:, 3], x[:, 4], x[:, 5]
N0, N1, N2 = x[:, 6], x[:, 7], x[:, 8]

X = (M0 - M) / M0  # 전환율(conversion)[cite: 11]

# 모멘트 합산[cite: 11]
mom0 = L0 + N0
mom1 = L1 + N1
mom2 = L2 + N2

# 5. 분자량 계산[cite: 11]
n = len(t)
Mn = np.zeros(n)
Mw = np.zeros(n)

# k=2부터 계산 (0으로 나누기 방지)[cite: 11]
for k in range(1, n):
    if mom0[k-1] != 0 and mom1[k-1] != 0:
        Mn[k] = mom1[k-1] / mom0[k-1]
        Mw[k] = mom2[k-1] / mom1[k-1]

PMn = 100.13 * Mn
PMw = 100.13 * Mw

# 6. 결과 시각화 (3x2 Subplot)[cite: 11]
plt.figure(figsize=(12, 15))

plt.subplot(3, 2, 1)
plt.plot(t, X), plt.xlabel('t(s)'), plt.ylabel('Conversion(X)')

plt.subplot(3, 2, 2)
plt.plot(t, I), plt.xlabel('t(s)'), plt.ylabel('Initiator(moles)')

plt.subplot(3, 2, 3)
plt.plot(t, M), plt.xlabel('t(s)'), plt.ylabel('Monomer(moles)')

plt.subplot(3, 2, 4)
plt.plot(t, Q), plt.xlabel('t(s)'), plt.ylabel('Q(kJ)')

plt.subplot(3, 2, 5)
plt.plot(t, PMn), plt.xlabel('t(s)'), plt.ylabel('Molecular weight PMn')

plt.subplot(3, 2, 6)
plt.plot(t, PMw), plt.xlabel('t(s)'), plt.ylabel('Molecular weight PMw')

plt.tight_layout()
plt.show()
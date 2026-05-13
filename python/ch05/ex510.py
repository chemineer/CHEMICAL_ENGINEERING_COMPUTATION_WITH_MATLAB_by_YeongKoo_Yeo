import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar # fsolve 대신 사용

# 상수 정의
di = 0.75e-2
z = 12
T = 418.15
P1 = 510
mu = 13.8e-6
M = 18
R = 8.314e-3

# 방정식 정의 (x가 P2)
def fP(x, G, P1, R, T, M, f, z, di):
    # P1^2 - x^2 = (G^2 * R * T / M) * (f * z / di + 2 * ln(P1/x))
    # 위 식을 0이 되도록 변환: x^2 + (G^2*R*T/M)*(f*z/di + 2*ln(P1/x)) - P1^2 = 0
    return x**2 + (G**2 * R * T / M) * (f * z / di + 2 * np.log(P1 / x)) - P1**2

G_range = np.arange(20, 50.1, 0.1)
P2_results = []

for G in G_range:
    Re = di * G / mu
    f = 1 / (0.79 * np.log(Re) - 0.64)**2
    
    # root_scalar(brentq) 사용: [0, P1] 구간 내에서 반드시 해를 찾음
    sol = root_scalar(fP, args=(G, P1, R, T, M, f, z, di), bracket=[0.001, P1])
    P2_results.append(sol.root)

# 결과 처리
P2_results = np.array(P2_results)
pressure_drop = P1 - P2_results

plt.figure(figsize=(8, 5))
plt.plot(G_range, pressure_drop)
plt.grid(True)
plt.xlabel('G (kg/sec/m^2)')
plt.ylabel('P_1 - P_2 (kPa)')
plt.show()
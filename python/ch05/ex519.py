import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# 상수 설정
Cv1 = 1.2e-2
Cv2 = Cv1
P0 = 1.014e5
P10 = 1.38e5
P30 = 1.08e5
P1 = P10
P3 = P30
z0_val = 3.05
rho = 1e3
A = 0.465
g = 9.8
P20 = P0 + rho * g * z0_val
cv = 0.2
R = 1.987
Mw = 29
V0 = 2.83
Vg0 = 1.415

k1 = np.sqrt(P10 / (P10 - P20))
k2 = k1 * Cv2 / Cv1
k3 = rho * g * z0_val / P10
k4 = V0 / Vg0
k5 = 0.735
k6 = R / (cv * Mw)

# 미분 방정식 정의
def cltank(z, t, k1, k2, k3, k4, k5, k6, P10, P1, P3):
    P1s = P1 / P10
    P3s = P3 / P10
    Vg = k4 - z
    Tg = (1.0 / Vg) ** k6
    Pg = k5 * Tg / Vg
    P2s = Pg + k3 * z
    F1s = k1 * np.sqrt(P1s - P2s)
    F2s = k2 * np.sqrt(P2s - P3s)
    return F1s - F2s

# 시간 범위 및 초기 조건
t = np.linspace(0, 0.4, 100)
z0 = 1.0

# 수치 적분 (ode45와 유사한 odeint 사용)
z = odeint(cltank, z0, t, args=(k1, k2, k3, k4, k5, k6, P10, P1, P3))

# 결과 계산
Vg = k4 - z.flatten()
Tg = (1.0 / Vg) ** k6
Pg = k5 * Tg / Vg
P1s = P1 / P10
P3s = P3 / P10
P2s = Pg + k3 * z.flatten()
F1s = k1 * np.sqrt(P1s - P2s)
F2s = k2 * np.sqrt(P2s - P3s)

# 그래프 출력
plt.plot(t, z, label='z')
plt.plot(t, F1s, '--', label='F_1')
plt.plot(t, F2s, '.-', label='F_2')
plt.plot(t, P2s, ':', label='P_2')
plt.legend()
plt.xlabel('t(dimensionless)')
plt.show()
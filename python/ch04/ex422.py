import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# --- 파라미터 설정 ---
T = 60 + 273.15  # Kelvin
R = 8.314        # J/(mol·K)
P1s = 83.25e3    # Pa
P2s = 37.97e3    # Pa
A12, A21 = 0.59, 1.42
B11, B22, B12 = -963e-6, -1523e-6, 52e-6  # m^3/mol

# --- 성분 계산 ---
x = np.linspace(0, 1, 101)  # x1 값 범위
gam1 = np.exp((1 - x)**2 * (A12 + 2 * (A21 - A12) * x))
gam2 = np.exp(x**2 * (A21 + 2 * (A12 - A21) * (1 - x)))
d12 = 2 * B12 - B11 - B22

# --- 비선형 방정식 시스템 정의 ---
def vle_system(z, x_val, g1, g2):
    """
    z[0] = y1, z[1] = P
    반환값: [방정식1, 방정식2]
    """
    y1, P = z
    # 기상 비이상성 보정 항 (Poynting 인자 포함 형태)
    phi1 = np.exp((B11 * (P - P1s) + P * d12 * (1 - y1)**2) / (R * T))
    phi2 = np.exp((B22 * (P - P2s) + P * d12 * y1**2) / (R * T))
    
    eq1 = x_val * g1 * P1s - y1 * P * phi1
    eq2 = (1 - x_val) * g2 * P2s - (1 - y1) * P * phi2
    return [eq1, eq2]

# --- 루프 계산 ---
y_res = []
P_res = []

for i in range(len(x)):
    # 초기 추정값 (y1=0.5, P=50000)
    z_sol = fsolve(vle_system, [0.5, 50000], args=(x[i], gam1[i], gam2[i]))
    y_res.append(z_sol[0])
    P_res.append(z_sol[1])

# --- 결과 시각화 ---
plt.plot(x, P_res, ':', label='Bubble point')
plt.plot(y_res, P_res, label='Dew point')
plt.xlabel('x1, y1')
plt.ylabel('P (Pa)')
plt.legend(loc='best')
plt.grid(True)
plt.show()
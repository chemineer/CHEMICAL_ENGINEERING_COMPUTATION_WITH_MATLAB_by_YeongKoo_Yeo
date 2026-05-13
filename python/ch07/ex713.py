import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 기본 데이터
x_data = np.array([0.0, 0.033, 0.072, 0.117, 0.171])
ye_data = np.array([0.0, 0.0396, 0.0829, 0.1127, 0.136])
V = 85; x0_init = 0.002; yp = 0.01

# (1) 커브 피팅
P3 = np.polyfit(x_data, ye_data, 3)
P4 = np.polyfit(x_data, ye_data, 4)
rmse3 = np.sqrt(np.mean((ye_data - np.polyval(P3, x_data))**2))
rmse4 = np.sqrt(np.mean((ye_data - np.polyval(P4, x_data))**2))
P = P4 if rmse4 < rmse3 else P3
print(f"Order of the fitting polynomial = {len(P)-1}")

xi = np.linspace(0, 0.18, 200); yi = np.polyval(P, xi)
ypr = yp / (1 - yp); x0r = x0_init / (1 - x0_init)

# (2) 조작선 플롯
L_list = [170, 150, 130]
plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.plot(xi, yi, label='Fitting curve'); plt.plot(x_data, ye_data, 'o', label='Data')
plt.xlabel('x'); plt.ylabel('y*'); plt.legend()

plt.subplot(2, 2, 2)
plt.plot(xi, yi, label='Equil. curve')
for L, style in zip(L_list, ['-', '--', ':']):
    r = L / V
    w = r * xi / (1 - xi) + ypr - r * x0r
    plt.plot(xi, w / (1 + w), style, label=f'L={L}')
plt.xlabel('x'); plt.ylabel('y'); plt.legend()

# (3) 이분법을 이용한 최소 유량(Lmin) 탐색
La, Lb = 80, 120; crit = abs(La - Lb)
def get_min_dist(L_val):
    r = L_val / V
    w = r * xi / (1 - xi) + ypr - r * x0r
    return np.min(w / (1 + w) - yi)

while crit > 1e-6:
    Lm = (La + Lb) / 2
    if get_min_dist(La) * get_min_dist(Lm) < 0:
        Lb = Lm
    else:
        La = Lm
    crit = abs(La - Lb)
Lmin = (La + Lb) / 2
print(f"Minimum liquid flow rate = {Lmin:.4f} kmol/h")

# (4) Lmin에서의 조작선 및 접점 계산
rm = Lmin / V
w_min = rm * xi / (1 - xi) + ypr - rm * x0r
plt.subplot(2, 2, 3)
plt.plot(xi, yi, label='Equil. curve'); plt.plot(xi, w_min / (1 + w_min), '--', label='L=Lmin')
plt.xlabel('x'); plt.ylabel('y'); plt.legend()

# 접점 x 찾기
def f_tangent(x):
    op = (rm * x / (1 - x) + ypr - rm * x0r) / (1 + rm * x / (1 - x) + ypr - rm * x0r)
    return op - np.polyval(P, x)

xmin = fsolve(f_tangent, 0.03)[0]
print(f"The operating curve becomes tangent at x = {xmin:.4f}")
plt.tight_layout(); plt.show()
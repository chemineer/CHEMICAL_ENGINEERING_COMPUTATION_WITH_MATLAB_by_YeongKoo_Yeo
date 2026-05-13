import numpy as np
import matplotlib.pyplot as plt

# 1. 데이터 정의 (온도 T 및 엔탈피 H)
T = np.array([283.15, 303.15, 323.15, 363.15, 393.15, 413.15]) # T(K)
H = np.array([2519.9, 2556.4, 2592.2, 2660.1, 2706.0, 2733.1]) # H(kJ/kg)

# 2. 2차 다항식 피팅 (2nd-order polynomial fitting)
# MATLAB: p = polyfit(T, H, 2)
p = np.polyfit(T, H, 2)

print("추정된 다항식 계수 (p2, p1, p0):")
print(p)

# 3. 보간 및 시각화 준비
Tv = np.arange(280, 415.1, 0.1)
# MATLAB: Hv = polyval(p, Tv)
Hv = np.polyval(p, Tv)

# 4. 특정 온도(350.15K)에서의 엔탈피 예측
# MATLAB: f = polyval(p, 350.15)
T_target = 350.15
f = np.polyval(p, T_target)

print(f"\n온도 {T_target} K에서의 예측 엔탈피:")
print(f"{f:.4f} kJ/kg")

# 5. 시각화
plt.figure(figsize=(8, 5))
plt.plot(Tv, Hv, label='2nd-order interpolation')
plt.plot(T, H, 'o', label='Steam table')
plt.xlabel('T(K)')
plt.ylabel('H(kJ/kg)')
plt.title('Enthalpy Interpolation for Saturated Steam')
plt.legend(loc='best')
plt.grid(True)
plt.show()
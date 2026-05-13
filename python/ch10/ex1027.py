import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin

# 1. 데이터 입력 (기질 농도 S, 반응 속도 r)
S = np.array([0.74, 1.19, 1.46, 1.62, 1.63, 5.27, 5.87, 6.26, 8.32, 10.10, 11.10,
              20.10, 21.73, 25.10, 27.78, 35.70])
r = np.array([0.10, 0.21, 0.21, 0.14, 0.34, 0.49, 0.36, 0.48, 0.82, 1.26, 0.55,
              3.34, 2.49, 2.01, 1.80, 1.68])

# 2. 초기 추정값 설정
rmax0 = 1.7
km0 = 11.0
x0 = [rmax0, km0]

# 3. 목적 함수(Cost Function) 정의: 오차 제곱합(SSE)
# x[0] = rmax, x[1] = km
def J(x):
    model = (x[0] * S) / (x[1] + S)
    return np.sum((r - model)**2)

# 4. 최적화 수행 (MATLAB의 fminsearch와 대응)
x_opt = fmin(J, x0)
rmax, km = x_opt

print(f"Maximum reaction rate (rm) = {rmax:.6g}")
print(f"Michaelis constant (km) = {km:.6g}")

# 5. 모델 비교 및 시각화
Si = np.linspace(np.floor(min(S)), np.ceil(max(S)), 1000)
ri = (rmax * Si) / (km + Si) # 최적화된 파라미터를 적용한 모델식

plt.figure(figsize=(8, 5))
plt.plot(Si, ri, label='Michaelis-Menten model') # 모델 선
plt.plot(S, r, 'o', label='Data')                # 실제 데이터 점
plt.xlabel('S(ng/ml)')
plt.ylabel('r(nmol/min)')
plt.legend(loc='best')
plt.grid(True)
plt.show()
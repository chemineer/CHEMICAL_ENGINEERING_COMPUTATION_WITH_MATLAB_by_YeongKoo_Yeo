import numpy as np
from scipy.optimize import least_squares, minimize_scalar

# 함수 정의: f(x) = 2*x^2*sin(x) + e^(-x)
f = lambda x: 2 * (x**2) * np.sin(x) + np.exp(-x)

# 1. lsqnonlin 방식 (비선형 최소 자승법)
# 초기값 x0 = -1, 하한 lb = -4, 상한 ub = 0
res_lsq = least_squares(f, x0=-1, bounds=(-4, 0))
x_lsq = res_lsq.x[0]

# 2. fminbnd 방식 (구간 내 최솟값 찾기)
# f(x)^2의 최솟값을 구함
f_sq = lambda x: (2 * (x**2) * np.sin(x) + np.exp(-x))**2
res_min = minimize_scalar(f_sq, bounds=(-4, 0), method='bounded')
x_min = res_min.x
fv_min = res_min.fun

# 결과 출력
print("--- lsqnonlin 변환 결과 ---")
print(f"x 값: {x_lsq:.6f}")

print("\n--- fminbnd 변환 결과 ---")
print(f"x 값: {x_min:.6f}")
print(f"f(x)^2의 최솟값: {fv_min:.6f}")
import numpy as np
from scipy.optimize import minimize_scalar

# 1. 목적 함수 정의
# f(x) = 2 * x^2 * sin(x) + e^(-x)
f = lambda x: 2 * (x**2) * np.sin(x) + np.exp(-x)

# 2. 최적화 수행 (구간 설정: lb = -4, ub = 0)
# 'bounded' 메서드를 사용하여 특정 구간 내의 최솟값을 찾습니다.
res = minimize_scalar(f, bounds=(-4, 0), method='bounded')

# 3. 결과 출력
print(f"최적의 x 값: {res.x}")
print(f"함수의 최솟값: {res.fun}")
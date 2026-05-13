import numpy as np
from scipy.integrate import dblquad

# 1. 2변수 함수 f(x, y) 정의
# scipy.dblquad는 함수의 인자 순서를 (y, x)로 받으므로 이에 맞춰 정의합니다.
fxy = lambda y, x: 3 * x**y + x - 1.2 * x**2 - 3 * y**2 + 25

# 2. 이중 적분 수행
# dblquad(func, x_low, x_high, y_low, y_high)
# x 범위: 0 ~ 7, y 범위: 0 ~ 5
area, error = dblquad(fxy, 0, 7, 0, 5)

# 3. 결과 출력
print(f"이중 적분 결과 (I2/Id): {area:.15f}")
print(f"추정 오차 (Estimated Error): {error:.15e}")
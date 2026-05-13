import numpy as np
from scipy.integrate import quad

# 1. 파라미터 정의
q = 0.2
r = 0.7
s = 4

# 2. 적분할 함수 hfun 정의
# x는 적분 변수입니다.
hfun = lambda x: 1 / ((x - q)**2 + 0.01) + 1 / ((x - r)**2 + 0.04) - s

# 3. 수치 적분 수행 (0부터 1까지)
# scipy.integrate.quad는 (적분값, 추정 오차)를 반환합니다.
result, error = quad(hfun, 0, 1)

# 4. 결과 출력
print(f"결과값 (Integral Result): {result:.15f}")
print(f"추정 오차 (Estimated Error): {error:.15e}")
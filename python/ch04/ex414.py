import numpy as np
from scipy.optimize import fsolve

# 1. 데이터 정의
x1 = np.array([0.0932, 0.1248, 0.1757, 0.2000, 0.2626, 0.3615, 0.4750, 0.5555, 0.6718])
Ge = np.array([-0.064, -0.086, -0.120, -0.133, -0.171, -0.212, -0.248, -0.252, -0.245])

# 2. 선형 회귀를 위한 행렬 X, Y 구성
# MATLAB: X = [x1.^2.*(1-x1); x1.*(1-x1).^2]'
Y = Ge.reshape(-1, 1)
col1 = (x1**2) * (1 - x1)
col2 = x1 * ((1 - x1)**2)
X = np.column_stack((col1, col2))

# 3. Margules 매개변수 계산 (최소자승법)
# MATLAB: A = inv(X'*X)*X'*Y
A = np.linalg.inv(X.T @ X) @ X.T @ Y
A21 = A[0][0]  # MATLAB의 A(1)
A12 = A[1][0]  # MATLAB의 A(2)

print(f"A12 = {A12:g}, A21 = {A21:g}")

# 4. gamma1 = gamma2가 되는 x1 찾기
# Margules 식 기반의 f(x) = exp(ln_gamma1) - exp(ln_gamma2)
def equation(x):
    ln_gamma1 = (1 - x)**2 * (A12 + 2 * (A21 - A12) * x)
    ln_gamma2 = x**2 * (A21 + 2 * (A12 - A21) * (1 - x))
    return np.exp(ln_gamma1) - np.exp(ln_gamma2)

x0 = 0.5
solution = fsolve(equation, x0)

print(f"찾은 x1 값: {solution[0]:.4f}")
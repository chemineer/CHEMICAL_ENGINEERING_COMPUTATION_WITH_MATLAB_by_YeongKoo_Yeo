import numpy as np
from scipy.optimize import minimize

# 1. 목적 함수 정의
# f(x) = (x1 - 1/2)^2 * (x1 + 1)^2 + 2 * (x2 + 1)^2 * (x2 - 1)^2
def objective(x):
    return (x[0] - 0.5)**2 * (x[0] + 1)**2 + 2 * (x[1] + 1)**2 * (x[1] - 1)**2

# 2. 제약 조건 설정 (A*x <= b)
# scipy에서 'ineq'는 fun(x) >= 0 을 의미하므로, b - A*x >= 0 으로 변형합니다.
A = np.array([[2, 4], 
              [-3, 1]])
b = np.array([7, 3])

def constraint(x):
    return b - np.dot(A, x)

cons = {'type': 'ineq', 'fun': constraint}

# 3. 최적화 수행
x0 = [0, 0]
res = minimize(objective, x0, method='SLSQP', constraints=cons)

# 4. 결과 출력
print(f"최적의 x 값: {res.x}")
print(f"함수의 최솟값: {res.fun}")
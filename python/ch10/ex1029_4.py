import numpy as np
from scipy.optimize import minimize

# 1. 행렬 및 벡터 정의
H = np.array([[2, -2], 
              [-2, 4]])
c = np.array([-4, 0])

# 부등식 제약 조건 Ax <= b
A = np.array([[2, 1], 
              [1, -4], 
              [-1, 0], 
              [0, -1]])
b = np.array([6, 0, 0, 0])

# 2. 목적 함수 정의 (1/2 * x.T * H * x + c.T * x)
def objective(x):
    return 0.5 * np.dot(x.T, np.dot(H, x)) + np.dot(c.T, x)

# 3. 제약 조건 설정 (b - Ax >= 0)
def constraint(x):
    return b - np.dot(A, x)

cons = {'type': 'ineq', 'fun': constraint}

# 4. 최적화 수행
x0 = np.zeros(2) # 초기값
res = minimize(objective, x0, method='SLSQP', constraints=cons)

# 5. 결과 출력
print(f"최적의 x 값: {res.x}")
print(f"함수의 최솟값: {res.fun}")
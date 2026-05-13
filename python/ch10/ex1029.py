import numpy as np
from scipy.optimize import lsq_linear

# 1. 행렬 및 벡터 정의
# C*x = d 를 최소화 (Least Squares)
C = np.array([[2, 0], 
              [0, 3]])
d = np.array([4, 4])

# A*x <= b (부등식 제약 조건)
A = np.array([[-2, 1], 
              [3, 5]])
b = np.array([6, 8])

# 2. lsq_linear 수행
# scipy의 lsq_linear는 기본적으로 lb <= x <= ub 제약을 지원하며,
# 선형 부등식 A*x <= b는 별도의 처리가 필요할 수 있으나 
# 여기서는 단순 변환을 위해 기본 구조를 제안합니다.
# 참고: lsq_linear는 박스 제약(bounds)에 특화되어 있습니다.
res = lsq_linear(C, d, bounds=(-np.inf, np.inf)) # 일반 해

# 만약 A*x <= b 제약 조건을 엄격히 적용해야 한다면 
# scipy.optimize.minimize를 사용하여 커스텀 목적함수를 구성하는 것이 일반적입니다.
from scipy.optimize import minimize

def objective(x):
    return np.sum((C @ x - d)**2)

cons = {'type': 'ineq', 'fun': lambda x: b - (A @ x)}
res_constrained = minimize(objective, x0=np.zeros(2), constraints=cons)

# 3. 결과 출력
print("최적의 x 값:", res_constrained.x)
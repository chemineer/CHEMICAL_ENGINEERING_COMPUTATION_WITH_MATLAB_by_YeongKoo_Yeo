import numpy as np
from scipy.optimize import minimize

# 1. 원본 함수 정의
# f(x) = 1/((x-0.3)^2+0.01) + 1/((x-0.9)^2+0.04) - 5
def f(x):
    return 1 / ((x - 0.3)**2 + 0.01) + 1 / ((x - 0.9)**2 + 0.04) - 5

# 2. Minimax를 위한 목적 함수 및 제약 조건 정의
# Minimax 문제는 'f(x) <= z'를 만족하는 최소의 z를 찾는 문제로 변환됩니다.
# x_aug[0]은 원래의 변수 x, x_aug[1]은 슬랙 변수 z입니다.

def objective(x_aug):
    return x_aug[1]  # z를 최소화

def constraint(x_aug):
    x = x_aug[0]
    z = x_aug[1]
    # z - f(x) >= 0 (즉, f(x) <= z)
    return z - f(x)

# 3. 최적화 수행
x0 = [1.0, f(1.0)]  # [초기 x값, 초기 z값(f(x0)로 설정)]
cons = {'type': 'ineq', 'fun': constraint}
res = minimize(objective, x0, method='SLSQP', constraints=cons)

# 4. 결과 출력
print(f"최적의 x 값: {res.x[0]:.6f}")
print(f"그때의 함수값 (Minimax): {res.x[1]:.6f}")
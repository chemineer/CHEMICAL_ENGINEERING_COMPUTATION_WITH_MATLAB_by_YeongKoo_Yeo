import numpy as np
from rosencopt import rosencopt

# 1. 목적 함수 및 그라디언트(delfun) 정의
# fun = @(x) 0.02*x(1)^2 + 1.2*x(2)^2 - 80
def fun(x):
    return 0.02 * (x[0]**2) + 1.2 * (x[1]**2) - 80

# delfun = @(x) [0.04*x(1); 2.4*x(2)]
def delfun(x):
    return np.array([0.04 * x[0], 2.4 * x[1]])

# 2. 문제 설정 (ex1013.m 내용 반영)
ne = 0                # 등식 제약 조건의 개수
m = 6                 # 총 제약 조건의 개수
crit = 1e-4           # 수렴 기준값
x0 = [4.0, 5.0]       # 초기값

# 제약 조건 행렬 A 및 벡터 b (Ax <= b 형태)
# 매트랩: A = [-1 0; -10 1; -1 0; 0 -1; 1 0; 0 1], b = [-3 -12 40 40 40 40]'
A = np.array([
    [-1.0,  0.0],
    [-10.0, 1.0],
    [-1.0,  0.0],
    [ 0.0, -1.0],
    [ 1.0,  0.0],
    [ 0.0,  1.0]
], dtype=float)

b = np.array([-3.0, -12.0, 40.0, 40.0, 40.0, 40.0], dtype=float)

# 3. rosencopt 함수 호출
# xopt: 최적해, fopt: 최적 목적 함수 값, iter_count: 반복 횟수
xopt, fopt, iter_count = rosencopt(fun, delfun, x0, A, b, ne, m, crit)

# 4. 결과 출력
print("--- Rosen's Gradient Projection Method 결과 ---")
print(f"최적해 (xopt): {xopt}")
print(f"최적값 (fopt): {fopt:.4f}")
print(f"반복 횟수 (iter): {iter_count}")
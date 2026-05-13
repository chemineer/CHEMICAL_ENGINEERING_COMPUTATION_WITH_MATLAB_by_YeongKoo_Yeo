import numpy as np
from zoutopt import zoutopt

# 1. 목적 함수 및 제약 조건 함수 정의 (funz)
# 매트랩의 funz 함수는 목적 함수(f)와 제약 조건(g)을 하나의 벡터로 반환합니다.
def funz(x):
    # f(x) = -(1.2*x1 + 3*x2)
    f = -(1.2 * x[0] + 3 * x[1])
    # g(x) = x1^2 + 6*x2^2 - 1 (g(x) <= 0 형태)
    g = np.array([x[0]**2 + 6 * x[1]**2 - 1])
    return f, g

# 2. 그라디언트 함수 정의
# delf: 목적 함수의 그라디언트
def delf(x):
    return np.array([-1.2, -3.0])

# delg: 제약 조건 함수의 그라디언트 (행렬 형태 또는 벡터 배열)
def delg(x):
    # g(x)의 x1, x2에 대한 편미분: [2*x1, 12*x2]
    return np.array([[2 * x[0], 12 * x[1]]])

# 3. 문제 설정 (ex1014.m 내용 반영)
x0 = np.array([1.0, 0.0])        # 초기값
xl = np.array([0.0, 0.0])        # 하한값 (Lower bounds)
xu = np.array([10.0, 10.0])      # 상한값 (Upper bounds)
nc = 0                 # 초기 활성 제약 조건 개수
ncs = 1                # 총 부등식 제약 조건 개수
crit = 1e-4            # 수렴 기준값
kmax = 1000            # 최대 반복 횟수

# 4. zoutopt 함수 호출
# xopt: 최적해, fopt: 최적 목적 함수 값, iter_count: 반복 횟수
xopt, fopt, iter_count = zoutopt(funz, delf, delg, x0, xl, xu, nc, ncs, crit, kmax)

# 5. 결과 출력
print("--- Zoutendijk's Method (ex1014) 결과 ---")
print(f"최적해 (xopt): {xopt}")
print(f"최적값 (fopt): {fopt:.4f}")
print(f"반복 횟수 (iter): {iter_count}")
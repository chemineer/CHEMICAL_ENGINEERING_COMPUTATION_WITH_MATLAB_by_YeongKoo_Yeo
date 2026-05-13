import numpy as np
from barnslp import barnslp

# 1. 문제 설정 (매트랩 ex1011.m 내용 반영)
# 제약 조건 행렬 A (이미 슬랙 변수가 포함된 형태)
A = np.array([[1, 1, 1, 1, 0], 
              [1, 2, 3, 0, 1]], dtype=float)

# 우변 상수 벡터 b
b = np.array([7, 12], dtype=float)

# 목적 함수 계수 c (Minimize f(x) = -2x1 - x2 - 4x3)
c = np.array([-2, -1, -4], dtype=float)

# 허용 오차(Tolerance)
tol = 1e-6

# 2. barnslp 함수 호출
# xopt: 최적해, fopt: 최적 목적 함수 값, basic: 기저 변수 인덱스
xopt, fopt, basic = barnslp(A, b, c, tol)

# 3. 결과 출력
print("--- Barnes 내점법(Interior Point Method) 결과 ---")
print(f"최적해 (xopt): {xopt}")
print(f"최소값 (fopt): {fopt:.4f}")
print(f"기저 변수 인덱스 (basic): {basic}")
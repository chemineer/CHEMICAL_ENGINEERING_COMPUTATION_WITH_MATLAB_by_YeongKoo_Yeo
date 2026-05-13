import numpy as np
from barnslp import barnslp

# 1. 문제 설정 (ex1012.m 내용 반영)
# A 행렬 초기화 (8행 13열)
A = np.zeros((8, 13))

# A(1:3, 1:5) 부분 설정
A[0:3, 0:5] = [
    [0.6, 0.5, 0.3, 0.4, 0.4],
    [0.2, 0.2, 0.3, 0.3, 0.1],
    [0.1, 0.2, 0.3, 0.2, 0.2]
]

# A(1:3, 6:8) 부분 설정 (3x3 단위행렬)
A[0:3, 5:8] = np.eye(3)

# A(4:8, 1:5) 부분 설정 (5x5 단위행렬)
A[3:8, 0:5] = np.eye(5)

# A(4:8, 9:13) 부분 설정 (5x5 단위행렬)
A[3:8, 8:13] = np.eye(5)

# 우변 상수 벡터 b (8x1)
b = 1e4 * np.array([17, 8.5, 7.5, 8, 10, 10, 10, 6])

# 목적 함수 계수 c (결정 변수 5개에 대한 계수)
c = np.array([-31.1, -20.7, -21.3, -23.2, -20.2])

# 허용 오차
tol = 1e-6

# 2. barnslp 함수 호출
# xopt: 최적해, fopt: 최적 목적 함수 값, basic: 기저 변수 인덱스
xopt, fopt, basic = barnslp(A, b, c, tol)

# 3. 결과 출력
print("--- Barnes 내점법(Interior Point Method) ex1012 결과 ---")
print(f"최적해 (xopt):\n{xopt}")
print(f"최소값 (fopt): {fopt:.4f}")
print(f"기저 변수 인덱스 (1-based): {basic}")
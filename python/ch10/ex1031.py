import numpy as np
from scipy.optimize import linprog

# 1. 목적 함수 계수 (c)
c = [-2, -3]

# 2. 부등식 제약 조건 (A*x <= b)
A = [
    [4, 10],
    [4, 4],
    [-1, 0],
    [0, -1]
]
b = [45, 23, 0, 0]

# 3. 정수 제약 조건 설정 (intcon = [1, 2])
# 모든 변수(x1, x2)가 정수여야 하므로 integrality에 1을 지정합니다.
# 1은 정수(Integer), 0은 연속 변수(Continuous)를 의미합니다.
integrality = [1, 1]

# 4. 정수 선형 계획법 수행
res = linprog(c, A_ub=A, b_ub=b, integrality=integrality, method='highs')

# 5. 결과 출력
if res.success:
    print(f"최적의 x 값: {res.x}")
    print(f"함수의 최솟값 (fv): {res.fun}")
else:
    print("최적해를 찾을 수 없습니다:", res.message)
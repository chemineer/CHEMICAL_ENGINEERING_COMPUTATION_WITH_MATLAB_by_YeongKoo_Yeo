import numpy as np
from scipy.optimize import linprog

# 1. 목적 함수 계수 (f)
# MATLAB: f = [-50 -98 -25 -43]'
c = [-50, -98, -25, -43]

# 2. 부등식 제약 조건 (A*x <= b)
A = [
    [6, 12, 3, 8],
    [4, 30, 2, 1],
    [-1, 0, 0, 0],
    [0, -1, 0, 0],
    [0, 0, -1, 0],
    [0, 0, 0, -1]
]
b = [1150, 750, 0, 0, 0, 0]

# 3. 선형 계획법 수행
# scipy.optimize.linprog는 기본적으로 목적 함수를 '최소화'합니다.
res = linprog(c, A_ub=A, b_ub=b, method='highs')

# 4. 결과 출력
if res.success:
    print(f"최적의 x 값: {res.x}")
    print(f"함수의 최솟값 (fv): {res.fun}")
else:
    print("최적해를 찾을 수 없습니다:", res.message)
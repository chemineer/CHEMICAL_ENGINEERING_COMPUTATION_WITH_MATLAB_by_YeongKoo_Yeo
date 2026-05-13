import numpy as np
from simplexlp import simplexlp

# 1. 문제 설정
A = [[2, 5], 
     [1, 1], 
     [3, 1]]
b = [20, 6, 9]
c = [2, 1]

# 2. 제약 조건 수정: 공백을 제거하여 '>>>'로 설정합니다.
# (simplexlp.py가 인덱스로 접근하므로 공백이 있으면 잘못된 조건으로 인식됩니다.)
constr = '>>>'

# 3. 함수 호출
result = simplexlp(A, b, c, constr)

# 4. 결과 출력
if result[0] is not None:
    xopt, fopt = result
    print("--- 심플렉스 알고리즘(Linear Programming) 결과 ---")
    print(f"최적해 (xopt): {xopt}")
    print(f"최소값 (fopt): {fopt:.4f}")
else:
    print("최적해를 찾을 수 없습니다.")
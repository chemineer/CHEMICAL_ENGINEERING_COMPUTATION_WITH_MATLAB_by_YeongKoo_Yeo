import numpy as np
from simplexlp import simplexlp

# 1. 문제 설정 (매트랩 ex1010.m 내용 반영)
# 목적 함수: Minimize f(x) = -x1 - x2
# 제약 조건:
# 0*x1 + 1*x2 <= 200
# 4*x1 + 3*x2 <= 1280
# 4*x1 + 1*x2 <= 960

A = [[0, 1], 
     [4, 3], 
     [4, 1]]
b = [200, 1280, 960]
c = [-1, -1]

# 2. 제약 조건 설정
# simplexlp.py의 내부 로직 처리를 위해 공백 없이 입력합니다.
constr = '<<<'

# 3. 심플렉스 함수 호출
result = simplexlp(A, b, c, constr)

# 4. 결과 출력
if result[0] is not None:
    xopt, fopt = result
    print("--- 심플렉스 알고리즘 (ex1010) 결과 ---")
    print(f"최적해 (xopt): {xopt}")
    print(f"최소값 (fopt): {fopt:.4f}")
    # 원래 문제가 최대화 문제였다면 (c가 음수이므로), 
    # 실제 최대값은 -fopt가 됩니다.
else:
    print("최적해를 찾을 수 없습니다.")
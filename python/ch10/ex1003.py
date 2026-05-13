import numpy as np
from brentopt import brentopt

# 1. 목적 함수 정의
# 매트랩의 @(x) exp(x) -3*x + 0.02/x - 0.00004/x 를 파이썬 함수로 변환
def fun(x):
    return np.exp(x) - 3*x + 0.02/x - 0.00004/x

# 2. 초기 매개변수 설정
x1 = 0.01
h = 0.2
crit = 1e-6

# 3. brentopt 함수 호출
xopt, fopt, nf = brentopt(fun, x1, h, crit)

# 4. 결과 출력
print("--- 최적화 결과 ---")
print(f"최적점 (xopt): {xopt:.8f}")
print(f"최솟값 (fopt): {fopt:.8f}")
print(f"함수 계산 횟수 (nf): {nf}")
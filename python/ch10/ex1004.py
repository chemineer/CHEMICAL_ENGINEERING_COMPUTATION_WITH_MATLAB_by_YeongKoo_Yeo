import numpy as np
from shubertopt import shubertopt

# 1. 목적 함수 정의
# 매트랩의 @(x) -sin(1.2*x)- sin(3.5*x) 를 파이썬 람다 함수로 변환
fun = lambda x: -np.sin(1.2 * x) - np.sin(3.5 * x)

# 2. 매개변수 설정
C = 8
a = -3
b = 8
crit = 1e-6
nfmax = 2000

# 3. shubertopt 함수 호출
xopt, fopt, nf = shubertopt(fun, a, b, C, crit, nfmax)

# 4. 결과 출력
print("--- Shubert-Piyavskii 최적화 결과 ---")
print(f"최적점 (xopt): {xopt:.8f}")
print(f"최대값 (fopt): {fopt:.8f}")
print(f"함수 계산 횟수 (nf): {nf}")
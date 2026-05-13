import numpy as np
from fbnopt import fbnopt  # fbnopt.py 모듈에서 함수 임포트

# 1. 초기 파라미터 설정
a = -8  # 구간 시작점
b = 8   # 구간 끝점
n = 20  # 반복 횟수

# 2. 목적 함수 정의
# MATLAB: fobj = @(x) 2*x^2*sin(1.5*x) - 4*x^2 + 3*x-1;
def fobj(x):
    return 2 * x**2 * np.sin(1.5 * x) - 4 * x**2 + 3 * x - 1

# 3. 피보나치 탐색 수행
x_opt, f_opt, fint = fbnopt(fobj, a, b, n)

# 4. 결과 출력
print(f"--- Fibonacci Search Result ---")
print(f"Optimal Point (x): {x_opt:.6f}")
print(f"Optimal Value (f): {f_opt:.6f}")
print(f"Final Interval Size: {fint:.6f}")
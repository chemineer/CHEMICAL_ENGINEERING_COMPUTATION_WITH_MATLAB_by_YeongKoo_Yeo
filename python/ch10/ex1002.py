import numpy as np
from gsopt import gsopt  # gsopt.py 모듈에서 함수 임포트

# 1. 초기 파라미터 설정
crit = 1e-6  # 허용 오차 (정지 조건)
h = 0.1      # 초기 보폭

# 2. 목적 함수 정의
# MATLAB: fobj = @(x) 2*x^2*sin(1.5*x) - 4*x^2 + 3*x-1;
def fobj(x):
    return 2 * x**2 * np.sin(1.5 * x) - 4 * x**2 + 3 * x - 1

# 3. 다양한 초기점(x1)에서 최적화 수행 및 결과 출력

# 케이스 1: x1 = -5
x1 = -5
x_opt, f_opt, n_eval = gsopt(fobj, x1, h, crit)
print(f"Case 1 (x1=-5): x = {x_opt:.6f}, f = {f_opt:.6f}, n = {n_eval}")

# 케이스 2: x1 = 1
x1 = 1
x_opt, f_opt, n_eval = gsopt(fobj, x1, h, crit)
print(f"Case 2 (x1= 1): x = {x_opt:.6f}, f = {f_opt:.6f}, n = {n_eval}")

# 케이스 3: x1 = 6
x1 = 6
x_opt, f_opt, n_eval = gsopt(fobj, x1, h, crit)
print(f"Case 3 (x1= 6): x = {x_opt:.6f}, f = {f_opt:.6f}, n = {n_eval}")
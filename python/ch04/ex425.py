import numpy as np
from scipy.optimize import fsolve
from beq import beq  # 모듈 import

# 1. 초기 데이터 설정
x110 = 0.0
x220 = 0.0
x210 = 1.0
x120 = 1.0
T0 = 100.0
beta0 = 0.8

z = np.array([0.2, 0.8])
P = 760.0
t0 = np.array([x110, x210, x120, x220, T0, beta0])

# 2. fsolve를 이용한 비선형 방정식 풀이
# MATLAB의 [x, fval] = fsolve(@beq, t0, [], z, P)에 대응
x, info, ier, msg = fsolve(beq, t0, args=(z, P), full_output=True)

if ier == 1:
    print("해를 찾았습니다.")
    print("결과값 (x11, x21, x12, x22, T, beta):", x)
else:
    print("해를 찾지 못했습니다:", msg)
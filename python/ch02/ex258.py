import numpy as np
from parapde import parapde

# 1. 초기 조건 및 경계 조건 설정
# f: 초기 온도 분포 (x^3 - 2*x^2 + 1.5*x)
# g0: 왼쪽 끝 경계 온도 (0)
# g1: 오른쪽 끝 경계 온도 (2)
f = lambda x: x**3 - 2 * x**2 + 1.5 * x
g0 = lambda t: np.zeros_like(t)
g1 = lambda t: np.full_like(t, 2)

# 2. 파라미터 설정
nx = 10
nt = 50
alpa = 0.8
tf = 0.1

# 3. parapde 함수 호출
# 수치적 안정성을 위해 r = alpa * dt / dx^2 값이 0.5 이하인지 확인이 필요합니다.
u = parapde(f, g0, g1, tf, nx, nt, alpa)
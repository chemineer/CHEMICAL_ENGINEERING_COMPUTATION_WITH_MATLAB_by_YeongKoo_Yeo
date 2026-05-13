import numpy as np
from hypbpde import hypbpde

# 1. 초기 조건 및 경계 조건 설정
# f1: 초기 변위 (x * (1 - x))
# f2: 초기 속도 (0)
# g0: 왼쪽 경계 조건 (0)
# g1: 오른쪽 경계 조건 (0)
f1 = lambda x: x * (1 - x)
f2 = lambda x: 0
g0 = lambda t: 0
g1 = lambda t: 0

# 2. 파라미터 설정
xspan = [0, 1]
tspan = [0, 1]
nx = 20
nt = 40
alpa = 1

# 3. hypbpde 함수 호출 (진동하는 현의 운동 계산 및 시각화)
# MATLAB의 [u r] = hypbpde(...) 구조를 유지
u, r = hypbpde(f1, f2, g0, g1, xspan, tspan, nx, nt, alpa)

# 4. 결과값 확인 (선택 사항)
print(f"Stability parameter (q): {r}")
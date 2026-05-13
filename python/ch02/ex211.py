import numpy as np
from scipy.optimize import fsolve

# 1. 함수 정의: x^3 - x*sin(x) + 1
f = lambda x: x**3 - x * np.sin(x) + 1

# 2. fsolve를 사용하여 x = -1 근처의 해를 찾음
# fsolve(함수, 초기 추측값)
x0 = -1
solution = fsolve(f, x0)

print(f"해(Root): {solution[0]:.6f}")
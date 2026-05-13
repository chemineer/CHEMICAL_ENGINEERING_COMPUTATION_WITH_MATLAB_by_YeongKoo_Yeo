import numpy as np
from trapzoid import trapzoid
from simps import simps

# 1. 함수 정의
f = lambda x: 1 / (1 + x**2)
a, b, n = 0, 1.5, 10

# 결과 실행
z_trap = trapzoid(f, a, b, n)
z_simp = simps(f, a, b, n)

print(f"Trapezoidal rule: {z_trap:.6f}")
print(f"Simpson 1/3 rule: {z_simp:.6f}")
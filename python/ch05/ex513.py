# ex513.py
import numpy as np
from scipy.optimize import fsolve
from pipsys import pipsys

# 1. 데이터 설정
L = 100 * np.array([6, 2, 1])
d = np.array([0.1524, 0.0625, 0.1016])
ep = 5e-6 * np.ones(3)
Z = np.array([8, 35, 50])
rho = 1e3
mu = 1e-3
g = 9.8
Q1 = 0.042

# 2. 초기값 설정
# x(0)=Ws, x(1)=v2, x(2)=v3, x(3)=f1, x(4)=f2, x(5)=f3
x0 = [10, 5, 5, 1e-3, 1e-3, 1e-3]

# 3. 비선형 방정식 시스템 풀이
# fsolve에 필요한 추가 인자들을 args로 전달
x = fsolve(pipsys, x0, args=(L, d, ep, Z, rho, mu, g, Q1))

Ws, v2, v3 = x[0], x[1], x[2]
Q2 = np.pi * v2 * (d[1]**2) / 4
Q3 = np.pi * v3 * (d[2]**2) / 4

# 4. 결과 출력
print(f"Power of the pump = {Ws:.4g} m of column of water")
print(f"Volumetric flow rate through pipe 1 = {Q1:.4g} m^3/sec")
print(f"Volumetric flow rate through pipe 2 = {Q2:.4g} m^3/sec")
print(f"Volumetric flow rate through pipe 3 = {Q3:.4g} m^3/sec")
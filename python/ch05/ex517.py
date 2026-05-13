import numpy as np
from scipy.optimize import fsolve
from pnetq import pnetq

# 매개변수 설정
L = 100 * np.array([1, 3, 12, 3, 12, 12, 3])
dP0 = -15e5
rho = 997.08
mu = 8.931e-4
D = 0.154
x0 = 0.1 * np.ones(7)

# fsolve 호출
# args 매개변수를 통해 pnetq에 추가 인자 전달
x = fsolve(pnetq, x0, args=(D, rho, mu, dP0, L))

print("계산된 유량 (x):", x)
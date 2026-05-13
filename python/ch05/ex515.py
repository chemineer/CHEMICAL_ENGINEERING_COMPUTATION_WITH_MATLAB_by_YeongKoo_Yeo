# ex515.py
import numpy as np
from scipy.optimize import fsolve
from dwfun import dwfun # 별도 파일에서 함수 임포트

# 1. 데이터 설정
W = 0.36
rho = 85.9
mu = 4.4e-4
rf = 5e-6
dP = 1e4
L = 4
K = 0.3
D0 = 0.1 # 초기값

# 2. 비선형 방정식 풀이
# fsolve(함수명, 초기값, args=(함수에 들어갈 나머지 인자들))
D = fsolve(dwfun, D0, args=(W, rho, mu, rf, dP, L, K))[0]

print(f"Calculated Internal Diameter (D): {D:.6f}")
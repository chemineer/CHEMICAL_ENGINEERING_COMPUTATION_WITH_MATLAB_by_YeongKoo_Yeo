# main.py
from scipy.optimize import fsolve
from nnDfun import nnDfun

# 데이터 설정
W = 6.67
rho = 961
dP = 1.5e4
rf = 5e-6
L = 10
k = 1.48
n = 0.64
K = 1.8

# 초기값 설정
D0 = 0.1

# fsolve 호출 (MATLAB의 fzero와 동일한 역할)
# fsolve는 배열을 반환하므로 [0]번 인덱스를 가져옵니다.
D = fsolve(nnDfun, D0, args=(W, rho, dP, rf, L, k, n, K))[0]

print(f"Calculated Pipe Diameter (D) = {D:.6f} m")
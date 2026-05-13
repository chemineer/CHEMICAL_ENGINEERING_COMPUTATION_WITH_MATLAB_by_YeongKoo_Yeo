# ex508.py
import numpy as np
from scipy.optimize import fsolve
from functools import partial
from vwfun import vwfun  # 별도 파일에서 함수 임포트

# 1. 초기값 및 상수 정의
T = 60
L = 1000
D = 7.981
rf = 0.00015
dz = 300
dP = -150
v0 = 10

# 2. 고정된 인자를 포함한 함수 생성
# fsolve는 첫 번째 인자로 변수(v)를 받는 함수만 처리 가능하므로 partial 사용
func = partial(vwfun, T=T, L=L, D=D, rf=rf, dz=dz, dP=dP)

# 3. 비선형 방정식 풀이
v = fsolve(func, v0)[0]

# 4. 결과 계산 (gpm 단위 변환)
# MATLAB: q = (7.481*60)*(pi*v.*(D/12).^2)/4
q = (7.481 * 60) * (np.pi * v * (D / 12)**2) / 4

print(f"Velocity (v): {v:.4f}")
print(f"Flow rate (q): {q:.4f} gpm")
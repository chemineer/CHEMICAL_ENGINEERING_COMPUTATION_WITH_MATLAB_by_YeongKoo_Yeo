import numpy as np
from scipy.optimize import fsolve
from dpf import dpf  # 모듈 import

# 1. 초기 데이터 설정
x10 = 0.2
x20 = 0.8
T0 = 100.0
y = np.array([0.2, 0.8])
P = 760.0

# 2. 초기 추정치 벡터
t0 = np.array([x10, x20, T0])

# 3. fsolve를 이용한 비선형 방정식 풀이
# MATLAB: x = fsolve(@dpf, t0, [], y, P)
x = fsolve(dpf, t0, args=(y, P))

print("계산된 결과 (x1, x2, T):")
print(x)
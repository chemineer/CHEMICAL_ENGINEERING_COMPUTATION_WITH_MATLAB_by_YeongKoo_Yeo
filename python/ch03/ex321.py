import numpy as np
from vpRM import vpRM

# 1. 물성치 및 파라미터 설정
# MATLAB: Tb = 159.22 + 273.15;
Tb = 159.22 + 273.15  # 끓는점 (K)

# MATLAB: nu = [2 1 1 5]; Dp = 0; GI = []; dB = [0.07545 -0.00227 0.11192 0.01653];
nu = np.array([2, 1, 1, 5])
Dp = 0
GI = []
dB = np.array([0.07545, -0.00227, 0.11192, 0.01653])

# 2. 증기압 계산 함수 호출 (100도 및 200도에서)
# MATLAB: T = 100 + 273.15; Pv = vpRM(T, Tb, nu, dB, GI, Dp);
T1 = 100 + 273.15
Pv1 = vpRM(T1, Tb, nu, dB, GI, Dp)
print(f"Vapor Pressure of n-Propylbenzene at {T1}K: {Pv1}")

# MATLAB: T = 200 + 273.15; Pv = vpRM(T, Tb, nu, dB, GI, Dp);
T2 = 200 + 273.15
Pv2 = vpRM(T2, Tb, nu, dB, GI, Dp)
print(f"Vapor Pressure of n-Propylbenzene at {T2}K: {Pv2}")
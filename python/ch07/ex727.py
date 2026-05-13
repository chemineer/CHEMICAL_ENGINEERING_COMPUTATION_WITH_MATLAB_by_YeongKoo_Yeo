import numpy as np

# 데이터 정의
t = np.array([4.4, 9.5, 16.3, 24.6, 34.7, 46.1, 59.0, 73.6, 89.4, 107.3]) # 시간 t(s)
V = 1e-3 * np.array([0.498, 1.0, 1.501, 2.0, 2.498, 3.002, 3.506, 4.004, 4.502, 5.009]) # 부피 V(m^3)

# 상수 정의
dP = 338e3      # 압력차 (Pa)
A = 0.0439      # 여과 면적 (m^2)
cs = 23.47      # 여과액 단위 부피당 케이크 질량 (kg/m^3)
mu = 8.937e-4   # 점도 (Pa·s)

# 선형 회귀 준비: t/V = (Kp/2) * V + B 형태
n = len(t)
Y = t / V
# MATLAB의 [V/2 ones(n,1)] 행렬 구성
C = np.column_stack((V / 2, np.ones(n)))

# 최소자승법 계산 (Normal Equation: X = (C^T * C)^-1 * C^T * Y)
X = np.linalg.inv(C.T @ C) @ C.T @ Y

Kp = X[0]
B = X[1]

# 물리적 파라미터 계산
alpha = Kp * (A**2) * dP / (mu * cs)
Rm = A * dP * B / mu

# 결과 출력
print(f"Alpha = {alpha:g} m/kg")
print(f"Rm = {Rm:g} m^(-1)")
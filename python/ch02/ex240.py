import numpy as np

# 1. 데이터 정의
t = np.array([4.4, 9.5, 16.3, 24.6, 34.7, 46.1, 59.0, 73.6, 89.4, 107.3])  # 시간 (sec)
V = 1e-3 * np.array([0.498, 1.0, 1.501, 2.0, 2.498, 3.002, 3.506, 4.004, 4.502, 5.009])  # 부피 (m^3)

A = 0.04        # 여과 면적
cs = 20         # 고형물 농도
vis = 8.937e-4  # 점도
dp = 3e5        # 압력차

# 2. 수치 미분 (dt/dV) 및 중점(Vm) 계산
dtV = np.diff(t) / np.diff(V)
n = len(V)
Vm = (V[:-1] + V[1:]) / 2

# 3. 선형 회귀 준비 (Y = C1*X1 + C2*X2)
# k1, k2 상수 정의
k1 = (vis * cs) / (A**2 * dp)
k2 = vis / (A * dp)

# 독립변수 행렬 X와 종속변수 벡터 Y 구성
# Y = [dt/dV]
# X = [k1*Vm, k2]
Y = dtV.reshape(-1, 1)
X = np.column_stack((k1 * Vm, k2 * np.ones(n-1)))

# 4. 정규 방정식을 이용한 회귀 분석 (MATLAB의 inv(X'*X)*X'*Y와 동일)
# 파이썬에서는 성능과 안정성을 위해 np.linalg.lstsq를 권장하지만, 원본 로직을 충실히 따름
C = np.linalg.inv(X.T @ X) @ X.T @ Y

# 5. 결과 출력
print(f"alpha = {C[0][0]:g}, Rm = {C[1][0]:g}")
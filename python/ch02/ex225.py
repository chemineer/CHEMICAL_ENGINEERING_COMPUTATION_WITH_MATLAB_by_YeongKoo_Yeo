import numpy as np
import matplotlib.pyplot as plt

# 1. 데이터 정의
T = np.array([25, 27, 30, 31, 35, 36, 37])
Pv = np.array([15.8, 26.43, 39.56, 94.91, 428.35, 861.71, 1851.24])

# 2. 선형 회귀를 위한 행렬 구성 (X * b = y 형태)
# MATLAB: X = [ones(1,length(T)); 1./T; log10(Pv)./T]'
logPv = np.log10(Pv)
X = np.column_stack([
    np.ones(len(T)), 
    1/T, 
    logPv/T
])

# 3. 정규 방정식을 이용한 계수 b 계산
# MATLAB: b = inv(X'*X)*X'*log10(Pv)'
b = np.linalg.inv(X.T @ X) @ X.T @ logPv

# 4. 안토안 계수 A, B, C 추출
A = b[0]
C = -b[2]
B = b[1] - A * C

print(f"추정된 계수:\n A = {A:.4f}\n B = {B:.4f}\n C = {C:.4f}")

# 5. 피팅 결과 시각화
Ti = np.arange(T[0], T[-1] + 0.1, 0.1)
Pvi = 10**(A + B / (Ti + C))

plt.figure(figsize=(8, 5))
plt.plot(Ti, Pvi, label='Fitting')
plt.plot(T, Pv, 'o', label='Data')
plt.xlabel('T(C)')
plt.ylabel('Pv(mmHg)')
plt.title('Vapor Pressure by Antoine Equation')
plt.legend()
plt.grid(True)
plt.show()
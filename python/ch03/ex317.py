import numpy as np
import matplotlib.pyplot as plt

# 1. 데이터 입력 (MATLAB과 동일)
t = np.array([0.001, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150])
Pv = np.array([65.4, 101.0, 159.3, 242.6, 357.9, 513.6, 718.7, 983.2, 1317.7, 1733.9, 
               2243.7, 2859.2, 3593.0, 4457.9, 5466.4, 6631.0])

# 2. 선형 회귀를 위한 행렬 구성
n = len(t)
Y = np.log(Pv).reshape(-1, 1)
# X = [1, 1/t, log(Pv)/t]
X = np.column_stack([np.ones(n), 1/t, np.log(Pv)/t])

# 3. K = inv(X'*X)*X'*Y 연산 (최소자승법)
K, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)

# 4. Antoine 상수 추출
A = K[0, 0]
C = -K[2, 0]
B = A * C - K[1, 0]

print(f"Antoine Constants: A={A:.4f}, B={B:.4f}, C={C:.4f}")

# 5. 그래프 생성을 위한 데이터 생성
ti = np.linspace(0, 150, 1500)
Pvi = np.exp(A - B / (ti + C))

# 6. 시각화
plt.figure(figsize=(8, 6))
plt.plot(ti, Pvi, label='Antoine eqn.')
plt.plot(t, Pv, 'o', label='Data')
plt.xlabel('t(deg.C)')
plt.ylabel('P_v(mmHg)')
plt.legend(loc='best')
plt.grid(True)
plt.show()
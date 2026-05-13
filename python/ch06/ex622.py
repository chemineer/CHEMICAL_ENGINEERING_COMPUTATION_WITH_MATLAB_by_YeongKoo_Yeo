import numpy as np
import matplotlib.pyplot as plt

# 1. 실험 데이터 입력
r = 1e-10 * np.array([71.0, 71.3, 41.6, 19.7, 42.0, 17.1, 71.8, 142.0, 284.0, 
                      47.0, 71.3, 117.0, 127.0, 131.0, 133.0, 41.8])
Pt = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0.5, 1, 5, 10, 15, 20, 1])
Ph = np.array([1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1, 1, 1, 1, 1, 1])
Pm = np.array([1, 4, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
Pb = np.array([0, 0, 1, 4, 1, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])

# 2. 선형 회귀를 위한 행렬 구성
# 반응 속도식 모델: r = (k * Ph * Pt) / (1 + Kb * Pb + Kt * Pt)
# 이를 선형 형태로 변환: (Ph * Pt / r) = (1/k) + (Kb/k) * Pb + (Kt/k) * Pt
n = len(r)
A = np.column_stack([np.ones(n), Pb, Pt])
b = (Ph * Pt) / r

# 3. 최소자승법 (Least Squares Method) 계산
# x = (A^T * A)^-1 * A^T * b
x = np.linalg.inv(A.T @ A) @ A.T @ b

# 4. 파라미터 추출
k = 1 / x[0]           # 반응 속도 상수
Kb = x[1] * k          # Pb에 대한 흡착 평형 상수
Kt = x[2] * k          # Pt에 대한 흡착 평형 상수

print(f"추정된 파라미터:")
print(f"k  = {k:.4e}")
print(f"Kb = {Kb:.4f}")
print(f"Kt = {Kt:.4f}")

# 5. 추정된 모델을 이용한 반응 속도 계산[cite: 20]
rc = (k * Ph * Pt) / (1 + Kb * Pb + Kt * Pt)
nc = np.arange(1, n + 1)

# 6. 결과 시각화[cite: 20]
plt.figure(figsize=(8, 6))
plt.plot(nc, r * 1e10, 'o', label='Data')
plt.plot(nc, rc * 1e10, '*', label='Estimated')
plt.xlabel('Run #')
plt.ylabel('$-10^{10}r$')
plt.title('Reaction Rate: Data vs Estimated')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
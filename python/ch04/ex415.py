import numpy as np
from unifgam import unifgam

# 데이터 정의
Nc = 2  # 성분의 수
k = 3   # 작용기의 수

# R: 각 작용기의 부피 벡터
R = np.array([0.9011, 0.6744, 1.2070])
# Q: 각 작용기의 표면적 벡터
Q = np.array([0.8480, 0.5400, 0.9360])

# nu: 작용기의 수 (행: 작용기 k, 열: 성분 i)
# MATLAB의 nu = [2 2; 1 5; 1 0]
nu = np.array([
    [2, 2],
    [1, 5],
    [1, 0]
])

# amn: 그룹 상호작용 매개변수 행렬
amn = np.array([
    [0, 0, 255.7],
    [0, 0, 255.7],
    [65.33, 65.33, 0]
])

T = 308.15      # 온도 (K)
x = np.array([0.4, 0.6])  # 각 성분의 몰 분율

# 활동도 계수 계산
gam = unifgam(k, R, Q, nu, amn, Nc, x, T)

# 결과 출력
print(f"Temperature: {T} K")
print(f"Mole Fractions: {x}")
print("-" * 30)
for i, g in enumerate(gam):
    print(f"Activity Coefficient gamma[{i+1}] = {g:.4f}")
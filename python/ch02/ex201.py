import numpy as np
from gausselm import gausselm
from GaussSeidel import GaussSeidel
from congrad import congrad

# --- 2. 메인 스크립트 실행 (ex201.m 내용) ---

# A 행렬 생성 (6x6)
# 대각 성분은 2, 인접 성분(i, i+1 및 i+1, i)은 -1인 삼중 대각 행렬
A = np.zeros((6, 6))
for i in range(6):
    A[i, i] = 2
for i in range(5):
    A[i, i+1] = -1
    A[i+1, i] = -1

# b 벡터 생성 및 초기값 설정
b = np.zeros(6)
b[5] = 4
rho = 1
x0 = np.zeros(6)

# 각 메서드별 계산 수행
xg = gausselm(A.copy(), b.copy()) # 원본 보존을 위해 copy 사용
xs = GaussSeidel(A, b, rho)
xc = congrad(A, x0, b)

# 결과 출력
print('x by Gauss elimination method:')
print(xg)

print('\nx by Gauss-Seidel method:')
print(xs)

print('\nx by conjugate gradient method:')
print(xc)